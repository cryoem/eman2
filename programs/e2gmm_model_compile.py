#!/usr/bin/env python
# Muyuan Chen 2023-10
from EMAN2 import *
import numpy as np
import protein_constant as e2pc
from Bio.PDB import *

floattype=np.float32
if "CUDA_VISIBLE_DEVICES" not in os.environ:
	# so we can decide which gpu to use with environmental variable
	os.environ["CUDA_VISIBLE_DEVICES"]='0' 
	
os.environ["TF_FORCE_GPU_ALLOW_GROWTH"]='true' 
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2' #### reduce log output
import tensorflow as tf
from e2gmm_model_refine import get_info, calc_bond, calc_angle, calc_dihedral_tf
	
def main():
	
	usage="""
	Compile stereochemical constraints from a pdb model that can be used for other parts of the GMM based refinement pipeline. 
	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--path", type=str,help="path for writing output", default=None)
	parser.add_argument("--model", type=str,help="load model from pdb file", default="")
	parser.add_argument("--writecif", action="store_true", default=False ,help="write cif instead of pdb.")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()
	logid=E2init(sys.argv,options.ppid)
	
	
	if options.model.endswith(".cif"):
		pdbpar = MMCIFParser( QUIET = True) 
	else:
		pdbpar = PDBParser( QUIET = True) 
		
	pdb = pdbpar.get_structure("model",options.model)
	
	######## remove Hydrogen if they are present in model
	residue=list(pdb.get_residues())
	nh=0
	for r in residue:
		d=list(r.child_dict)
		for a in d:
			if a[0]=='H':
				r.detach_child(a)
				nh+1
	if nh>0:
		print(f"Removing {nh} H atoms")
	atoms=list(pdb.get_atoms())
	
	##########################
	atom_pos=np.array([a.get_coord() for a in atoms])
	atom_res=np.array([a.get_parent().get_id()[1] for a in atoms])
	atom_chn=np.array([a.get_parent().get_parent().get_id() for a in atoms])
	c,cid=np.unique(atom_chn, return_inverse=True)
	
	if options.path==None: 
		path=options.path=num_path_new("gmm_model_")
	else:
		path=options.path
		
	##########################
	options.cmd=' '.join(sys.argv)
	fm=f"{options.path}/0_compile_params.json"
	js=js_open_dict(fm)
	js.update(vars(options))
	js.close()	
	
	if options.model.endswith(".cif") or options.writecif:
		io=MMCIFIO()
		io.set_structure(pdb)
		io.save(f"{path}/model_input.cif")
	else:
		io=PDBIO()
		io.set_structure(pdb)
		io.save(f"{path}/model_input.pdb")
		
	##########################
	print("Calculating van der Waals radius...")
	atomtype=['_'.join([a.parent.resname, a.id]) for a in atoms]
	vdw=np.array([e2pc.get_vdw_radius(a) for a in atomtype])
	np.savetxt(f"{path}/model_vdwr.txt", vdw)
	
	
	##########################
	print("Compiling bonds...")
	bonds0=[] # intra-residue bonds
	bonds1=[] # peptide bonds
	bonds2=[] # S-S bonds
	bonds3=[] # bonds between DNA/RNA base paris 
	sgid=np.array([i for i,a in enumerate(atoms) if a.get_id()=='SG'])

	ncbond=[1.336, 0.012] ## peptide bonds length and std
	ssbond=[2.033, 0.04] ## S-S bonds length and std
	dnabond=[1.601, 0.012] ## DNA/RNA bonds length and std
	for ai,at in enumerate(atoms):
		idx_res=(atom_chn==atom_chn[ai])*(atom_res==atom_res[ai])
		idx_res=np.where(idx_res)[0]
		resname=at.get_parent().get_resname()
		atname=at.get_id()
		for i2 in idx_res:
			a2=atoms[i2]
			a2name=a2.get_id()
			ky0='_'.join([resname, atname, a2name])
			if ky0 in e2pc.bond_length_std:
				sc=e2pc.bond_length_std[ky0]
				bonds0.append([ai, i2, sc[0], sc[1]])
				
		if atname=='N':
			idx_res=(atom_chn==atom_chn[ai])*(atom_res==atom_res[ai]-1)
			idx_res=np.where(idx_res)[0]
			idx_res=[i for i in idx_res if atoms[i].get_id()=='C']
			if len(idx_res)==0: continue
			i2=idx_res[0]
			bonds1.append([ai, i2, ncbond[0], ncbond[1]])
			
		if atname=='SG':
			dp=atom_pos[sgid]-atom_pos[ai]
			dp=np.linalg.norm(dp, axis=1)
			dp=abs(dp-2.034)
			i2=np.argmin(dp)
			dp=dp[i2]
			i2=sgid[i2]
			if dp<.1 and i2>ai:
				bonds2.append([ai, i2, ssbond[0], ssbond[1]])
				
		if atname=="O3'":
			idx_res=(atom_chn==atom_chn[ai])*(atom_res==atom_res[ai]+1)
			idx_res=np.where(idx_res)[0]
			idx_res=[i for i in idx_res if atoms[i].get_id()=='P']
			if len(idx_res)==0: continue
			i2=idx_res[0]
			bonds3.append([ai, i2, dnabond[0], dnabond[1]])
			
	bonds0=np.array(bonds0).reshape(-1,4)
	bonds1=np.array(bonds1).reshape(-1,4)
	bonds2=np.array(bonds2).reshape(-1,4)
	bonds3=np.array(bonds3).reshape(-1,4)
	print("  {} intra-residue bonds, {} peptide bonds, {} DNA/RNA bonds, and {} S-S bonds.".format(len(bonds0), len(bonds1), len(bonds3), len(bonds2)))

	bonds=np.concatenate([bonds0, bonds1, bonds2, bonds3], axis=0)
	print("  {} bonds total.".format(len(bonds)))
	
	bond_len=calc_bond(atom_pos[None,:,:], bonds[:,:2].astype(int))
	ii=np.where(abs(bond_len[0])>15)[0]
	if len(ii)>0:
		print("  ignore bad bonds:")
	bd=bonds[ii][:,:2].astype(int)
	for i,b in enumerate(bd):
		bx=[get_info(atoms[x], include_id=True) for x in b]
		print("    {} - {} : length = {:.1f} A".format(bx[0], bx[1], float(bond_len[0,ii[i]])))
	bonds=np.delete(bonds, ii, axis=0)
	
	bond_len=calc_bond(atom_pos[None,:,:], bonds[:,:2].astype(int))
	bond_df=(bond_len-bonds[:,2])/bonds[:,3]
	bond_df=abs(bond_df[0].numpy())
	print("  Average bond length deviation: {:.2f} std. {} outliers beyond 5 std.".format(np.mean(bond_df), np.sum(bond_df>5)))
	np.savetxt(f"{path}/model_bond.txt", bonds)
		
		
	##########################
	print("Calculating bond angles...")
	angs=[]
	skip=[]
	for ai,at in enumerate(atoms):
		nb0=bonds[bonds[:,0]==ai,1]
		nb1=bonds[bonds[:,1]==ai,0]
		nbs=np.append(nb0, nb1).astype(int)
		if len(nbs)<2: continue
		resname=atoms[ai].get_parent().get_resname()
		pairs=[(i0,i1) for i0 in nbs for i1 in nbs if i0!=i1]
		for i0,i1 in pairs:
			ky0='_'.join([atoms[i].get_id() for i in [i0,ai,i1]])
			ky='_'.join([resname, ky0])
			if ky=="LEU_C_N_CA": 
				print(ky, [atoms[i].parent.id[1] for i in [i0,ai,i1]])
			if ky in e2pc.bond_angle_std:
				sc=e2pc.bond_angle_std[ky]
				angs.append([i0,ai,i1, sc[0], sc[1]])
			else:
				ky1='_'.join([atoms[i].get_id() for i in [i1,ai,i0]])
				kyr='_'.join([resname, ky1])
				if not kyr in e2pc.bond_angle_std:
					skip.append(ky)
		
	angs=np.array(angs)
	if len(skip)>0:
		print("Skipping following angles:")
		for k in skip:
			print(k)
	
	ang_val=calc_angle(atom_pos[None,:,:], angs[:,:3].astype(int))
	ang_df=(ang_val-angs[:,3])/angs[:,4]
	ang_df=abs(ang_df[0].numpy())
	print("  {} bonds angles total.".format(len(angs)))
	print("  Average bond angle deviation: {:.2f} std. {} outliers beyond 5 std.".format(np.mean(ang_df), np.sum(ang_df>5)))
	
	np.savetxt(f"{path}/model_angle.txt", angs)
	
	
	##########################
	print("Compiling Ramachandran angles...")
	idx_dih=[]
	atomtype=np.array([a.get_id() for a in atoms])
	thr=1.6
	for ci in np.unique(atom_chn):
		pcid=np.where(atom_chn==ci)[0]
		ares=atom_res[pcid]
		for ri in np.unique(atom_res):
			i0=pcid[ares==ri]
			pr=np.hstack([atom_pos[i0],i0[:,None]])
			
			i1=pcid[(ares==ri-1) + (ares==ri+1)]
			pr1=np.hstack([atom_pos[i1],i1[:,None]])
			if len(pr1)==0: continue
			
			pa_n=pr[atomtype[i0]=='N']
			pa_c=pr[atomtype[i0]=='C']
			pa_ca=pr[atomtype[i0]=='CA']
			if len(pa_n)!=1 or len(pa_c)!=1 or len(pa_ca)!=1: continue
			
			dst_n=np.linalg.norm((pr1-pa_n)[:,:3], axis=1)
			dst_c=np.linalg.norm((pr1-pa_c)[:,:3], axis=1)
			if np.min(dst_n)>thr or np.min(dst_c)>thr: continue
			
			pa_nb0=pr1[np.argmin(dst_n)]
			pa_nb1=pr1[np.argmin(dst_c)]
			
			phi=np.vstack([pa_nb0, pa_n, pa_ca, pa_c])
			psi=np.vstack([pa_n, pa_ca, pa_c, pa_nb1])
			
			idx=np.append(phi[:,3], psi[:,3])
			idx_dih.append(idx)
			
	idx_dih=np.array(idx_dih, dtype=int)
	print("  {} ramachandran angles total.".format(len(idx_dih)))
	
	np.savetxt(f"{path}/model_rama_angle.txt", idx_dih)
	
	##########################
	print("Compiling dihedral angles...")
	
	a=np.diff(atom_res)
	a=(a!=0).astype(int)
	a=np.append(0,a)
	atom_residx=np.cumsum(a)
	ca_res=np.array([a.get_parent().get_resname() for a in atoms])
	a,i=np.unique(atom_residx, return_index=True)
	ca_res=ca_res[i]
	
	dihs=[]
	dihs_type=[]
	aname_last=[]
	for ii in np.unique(atom_residx):
		res=ca_res[ii]
		if res not in e2pc.residue_dihedral: continue
		dh=e2pc.residue_dihedral[res]

		ia=np.where(atom_residx==ii)[0]
		aname={atoms[i].get_id():i for i in ia}
			
		for h in dh:
			chid=[aname[a] for a in h[1:5] if a in aname]
			if len(chid)<4:
				continue
			dihs.append(chid)
			dihs_type.append(h[0])
		
		if len(aname_last)>0:
			
			if ('CA' in aname_last) and  ('CA' in aname):
				r0=atoms[aname_last['CA']].get_parent().get_id()[1]
				r1=atoms[aname['CA']].get_parent().get_id()[1]
				c0=atoms[aname_last['CA']].get_parent().get_parent().get_id()
				c1=atoms[aname['CA']].get_parent().get_parent().get_id()
				if r0+1==r1 and c0==c1:
					try:
						dh0=[aname_last['CA'], aname_last['C'], aname['N'], aname['CA']]
						dh1=[aname_last['O'], aname_last['C'], aname_last['CA'], aname['N']]
					except:
						print(f"missing atoms in {ii} {res}")
						continue
					dihs.append(dh0)
					dihs.append(dh1)
					dihs_type.extend(['peptide','backbone'])
					
		aname_last=aname
		
	dihs=np.array(dihs, dtype=int)
	dihs_type=np.array(dihs_type)
	print("  {} planar dihedral angles.".format(len(dihs)))
	
	pid=[i for i,d in enumerate(dihs) if dihs_type[i]=="backbone"]
	dihs_backbone=dihs[pid]
	print("    backbone:  {}".format(len(dihs_backbone)))

	pid=[i for i,d in enumerate(dihs) if dihs_type[i]=="peptide"]
	dihs_peptide=dihs[pid]
	print("    peptide:   {}".format(len(dihs_peptide)))
	
	pid=[i for i,d in enumerate(dihs) if dihs_type[i].startswith("C")]
	dihs_sidechain=dihs[pid]
	print("    sidechain: {}".format(len(dihs_sidechain)))
	
	dihs_plane=np.vstack([dihs_backbone, dihs_sidechain])
	pt=tf.gather(atom_pos[None,:], dihs_plane[:,:4], axis=1)
	rot=calc_dihedral_tf(pt)[0].numpy()
	rot=np.sin(np.deg2rad(rot))
	rot=np.degrees(np.arcsin(abs(rot)))
	print("  Average backbone/sidechain planar dihedral angle deviation: {:.2f} degrees".format(np.mean(rot)))
	print("  {} angle outliers beyond 10 degrees".format(np.sum(rot>10)))
	# ir=np.where(rot>10)[0]
	# for i in ir:
	# 	print(get_info(atoms[dihs_plane[i,2]]), rot[i])
	
	if len(dihs_peptide)>0:
		pt=tf.gather(atom_pos[None,:], dihs_peptide[:,:4], axis=1)
		rot=calc_dihedral_tf(pt)[0].numpy()
		rot=np.sin(np.deg2rad(rot))
		rot=np.degrees(np.arcsin(abs(rot)))
		print("  Average peptide dihedral angle deviation: {:.2f} degrees".format(np.mean(rot)))
		print("  {} angle outliers beyond 30 degrees".format(np.sum(rot>30)))
		#ii=np.where(rot>30)[0]
		#for i in ii:
			#d=dihs_peptide[i,:4]
			#print(i, rot[i], [atoms[j].get_parent().get_id()[1] for j in d])
	
	np.savetxt(f"{path}/model_dih_plane.txt", dihs_plane.astype(int))
	np.savetxt(f"{path}/model_dih_piptide.txt", dihs_peptide.astype(int))
	
	chis=[f"chi{i}" for i in range(1,5)]
	dihs_chi_id=list([i for i,d in enumerate(dihs) if dihs_type[i] in chis])
	dihs_chi=dihs[dihs_chi_id]
	dihs_chi_res=[atoms[i].get_parent().get_resname() for i in dihs_chi[:,1]]
	# for i in dihs_chi[:,1]:
	# 	print(get_info(atoms[i]))
	dihs_chi_res=np.array([e2pc.restype_3_to_index[i] for i in dihs_chi_res])
	chi_id=np.array([int(i[3]) for i in dihs_type[dihs_chi_id]])
	print("{} sidechain rotamer chi angles".format(len(chi_id)))

	tosave=np.hstack([dihs_chi, dihs_chi_res[:,None], chi_id[:,None]])
	np.savetxt(f"{path}/model_dih_chi.txt", tosave.astype(int))
	
	#### RNA backbone...
	residues=list(pdb.get_residues())
	rid=[r.resname for r in residues]
	if np.any([a in rid for a in "AUGC"]):
		print("RNA detected")

		# bkb=["P", "O5'", "C5'", "C4'", "C3'", "O3'"]
		idx_rna_dih=[]
		for ir, res in enumerate(residues):
			# print(ir, res)
			if ir==0: continue
			dic0=residues[ir-1].child_dict
			dic1=res.child_dict
			try:
				b0=[dic0[b] for b in ["C5'", "C4'", "C3'", "O3'"]]
				b1=[dic1[b] for b in ["P", "O5'", "C5'", "C4'", "C3'", "O3'"]]
			except:
				print("missing atoms at",res)
				continue
			
			bk=b0+b1
			bk=[atoms.index(b) for b in bk]
			
			dih=[bk[i:i+4] for i in range(7)]
			
			idx_rna_dih.append(dih)
			
		idx_rna_dih=np.array(idx_rna_dih)
		# print(idx_rna_dih.shape)
		idx_rna_dih=idx_rna_dih.reshape((idx_rna_dih.shape[0], -1))
		print("{} RNA backbone angles".format(len(idx_rna_dih)))

		# tosave=np.hstack([dihs_rna, dihs_rna_res[:,None], rna_id[:,None]])
		np.savetxt(f"{path}/model_dih_rna.txt", idx_rna_dih.astype(int))
	
	print(f"Done. Parameters written in folder {path}")
	E2end(logid)
	
	
if __name__ == '__main__':
	main()
	
