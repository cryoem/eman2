
def generate_das_without_psi(delta, output_filename=None, method='S', symmetry="c1"):
    from utilities import even_angles
    from utilities import write_text_row    
    a = even_angles(method=method, delta=delta, theta2=180.0, symmetry=symmetry)
    for i in a:
        i[2] = 0.0
        while len(i) < 5:
            i.append(0.0)
    if output_filename != None and output_filename != "":
        write_text_row(a, output_filename)
    return a
       
       
def generate_das_with_psi(delta, output_filename=None, method='S', symmetry="c1"):
    from utilities import even_angles
    #from utilities import write_text_row    
    a = even_angles(method=method, delta=delta, theta2=180.0, symmetry=symmetry)
    steps_count = int(round(360.0/delta))
    step = 360.0 / steps_count
    a2 = []
    for i in a:
        for j in xrange(steps_count):
            new_params = i[:]
            new_params[2] = j * step
            while len(new_params) < 5:
                new_params.append(0.0)
            a2.append(new_params)
    if output_filename != None and output_filename != "":
        #write_text_row(a2, output_filename)
        outf = open(output_filename, "w")
        for i in xrange(len(a2)):
            for j in xrange(len(a2[i])):
                outf.write(" %.5g"%a2[i][j])
            outf.write("\n")
        outf.close()
    return a2
        
        
def generate_das_full(delta, output_filename=None, method='S', symmetry="c1"):
    from utilities import even_angles
    #from utilities import write_text_row    
    a = even_angles(method=method, delta=delta, theta2=180.0, symmetry=symmetry)
    steps_count = int(round(360.0/delta))
    step = 360.0 / steps_count
    a2 = []
    for i in a:
        for j in xrange(steps_count):
            new_params = i[:]
            new_params[2] = j * step
            while len(new_params) < 5:
                new_params.append(0.0)
            for sx in [-2, -1, 0, 1, 2]:
                for sy in [-2, -1, 0, 1, 2]:
                    new_params2 = new_params[:]
                    new_params2[3] = sx
                    new_params2[4] = sy
                    a2.append(new_params2)
    if output_filename != None and output_filename != "":
        #write_text_row(a2, output_filename)
        outf = open(output_filename, "w")
        for i in xrange(len(a2)):
            for j in xrange(len(a2[i])):
                outf.write(" %.5g"%a2[i][j])
            outf.write("\n")
        outf.close()
    return a2
           

def generate_projections_directions(count, half_sphere=False, output_filename=None, method='S', symmetry="c1"):  # method=P|R|S (R - random)
    """
    Generates set of projections directions with given method and size (psi and shifts are set to 0)
    Returns: list of projections directions, delta
    """
    from utilities import even_angles
    from random import random, uniform
    from utilities import write_text_row
    from math import sin, pi
    
    max_theta = 180.0
    if half_sphere:
        max_theta = 90.0
    
    if method == 'R':
        a1 = []
        while len(a1) < count:
            theta = uniform(0.0, max_theta)
            if random() < abs(sin(theta)):
                a1.append( [uniform(0.0, 360.0), theta/pi*180.0] )  
    
    if method == 'P' or method == 'S':
        delta1 = 10.0
        delta2 = 5.0
        a1 = even_angles(method=method, delta=delta1, theta2=max_theta, symmetry=symmetry)
        a2 = even_angles(method=method, delta=delta2, theta2=max_theta, symmetry=symmetry)
        while len(a1) > count:
            delta2 = delta1
            a2 = a1
            delta1 *= 2.0
            a1 = even_angles(method=method, delta=delta1, theta2=max_theta, symmetry=symmetry)
            #print "1", delta1, delta2
        while len(a2) < count:
            delta1 = delta2
            a1 = a2
            delta2 /= 2.0
            a2 = even_angles(method=method, delta=delta2, theta2=max_theta, symmetry=symmetry)
            #print "2", delta1, delta2
        while len(a1) < count:
            delta3 = (delta1 + delta2) / 2.0
            a3 = even_angles(method=method, delta=delta3, theta2=max_theta, symmetry=symmetry)
            if len(a3) > count:
                if delta2 == delta3: break
                delta2 = delta3
                a2 = a3
            else:
                if delta1 == delta3: break
                delta1 = delta3
                a1 = a3
            #print "3", delta1, delta2
        #print "delta=", delta1
    
    for pos in a1:
        for i in range(2,5):
            if len(pos) <= i:
                pos.append(0.0)
            else:
                pos[i] = 0.0
    
    if output_filename != None and output_filename != "":
        write_text_row(a1, output_filename)
    
    return a1, delta1


def mirror_projection_direction(pAngles):
    angles = pAngles[:]
    angles[0] = angles[0] + 180.0
    angles[1] = 180.0 - angles[1]
    if angles[0] >= 360.0: angles[0] -= 360.0
    angles[2] = -angles[2]
    if angles[2] <    0.0: angles[2] += 360.0    
    return angles


def prepare_projections_without_psi(volume, count, output_filename=None, tanl_filter_cutoff=0.35, tanl_filter_fall_off=0.2, symmetry="c1"):
    """
    Prepare projections
    Returns: list of projections
    """
    from EMAN2 import EMData
    from filter import filt_tanl
    from projection import project
    from random import random
    import types
    
    if type(volume) is types.StringType:
        v = EMData()
        v.read_image(volume)
        volume = v  
    
    proj_directions = generate_projections_directions(count, True, symmetry=symmetry)[0]
    
    nx = volume.get_xsize()
    projections = [] 
    for iProj in xrange(len(proj_directions)):
        if random() < 0.5:
            proj_directions[iProj] = mirror_projection_direction(proj_directions[iProj])
        projections.append( project(volume, proj_directions[iProj], nx // 2 - 1) )
        projections[iProj] = filt_tanl(projections[iProj], tanl_filter_cutoff, tanl_filter_fall_off)
        
    if output_filename != None and output_filename != "":
        for iProj in xrange(len(projections)):
            projections[iProj].write_image(output_filename, iProj)
        
    return projections


def prepare_projections_with_psi(volume, count, output_filename=None, tanl_filter_cutoff=0.35, tanl_filter_fall_off=0.2, symmetry="c1"):
    """
    Prepare projections
    Returns: list of projections
    """
    from EMAN2 import EMData
    from filter import filt_tanl
    from projection import project
    from random import random
    import types
    
    if type(volume) is types.StringType:
        v = EMData()
        v.read_image(volume)
        volume = v  
    
    proj_directions = generate_projections_directions(count, True, symmetry=symmetry)[0]
    
    nx = volume.get_xsize()
    projections = [] 
    for iProj in xrange(len(proj_directions)):
        proj_directions[iProj][2] = random() * 360.0
        if random() < 0.5:
            proj_directions[iProj] = mirror_projection_direction(proj_directions[iProj])
        projections.append( project(volume, proj_directions[iProj], nx // 2 - 1) )
        projections[iProj] = filt_tanl(projections[iProj], tanl_filter_cutoff, tanl_filter_fall_off)
        
    if output_filename != None and output_filename != "":
        for iProj in xrange(len(projections)):
            projections[iProj].write_image(output_filename, iProj)
        
    return projections


def prepare_projections_full(volume, count, output_filename=None, tanl_filter_cutoff=0.35, tanl_filter_fall_off=0.2, symmetry="c1"):
    """
    Prepare projections
    Returns: list of projections
    """
    from EMAN2 import EMData
    from filter import filt_tanl
    from projection import project
    from random import random, shuffle
    from utilities import even_angles
    import types
    
    if type(volume) is types.StringType:
        v = EMData()
        v.read_image(volume)
        volume = v  
    
    proj_directions = even_angles(method='S', delta=1.0, theta2=180.0, symmetry=symmetry)
    shuffle(proj_directions)
    proj_directions = proj_directions[:count]
    
    nx = volume.get_xsize()
    projections = [] 
    for iProj in xrange(len(proj_directions)):
        while len(proj_directions[iProj]) < 5:
            proj_directions[iProj].append(0)
        proj_directions[iProj][2] = random() * 360.0
        proj_directions[iProj][3] = random() * 4.0 - 2.0
        proj_directions[iProj][4] = random() * 4.0 - 2.0
        projections.append( project(volume, proj_directions[iProj], nx // 2 - 1) )
        projections[iProj] = filt_tanl(projections[iProj], tanl_filter_cutoff, tanl_filter_fall_off)
        
    if output_filename != None and output_filename != "":
        for iProj in xrange(len(projections)):
            projections[iProj].write_image(output_filename, iProj)
        
    return projections


#def prepare_projections_full(volume, count, output_filename=None, tanl_filter_cutoff=0.35, tanl_filter_fall_off=0.2, symmetry="c1"):
#    """
#    Prepare projections
#    Returns: list of projections
#    """
#    from EMAN2 import EMData
#    from filter import filt_tanl
#    from projection import project
#    from random import random, randrange
#    import types
#    
#    if type(volume) is types.StringType:
#        v = EMData()
#        v.read_image(volume)
#        volume = v  
#    
#    proj_directions = generate_projections_directions(count, True, symmetry=symmetry)[0]
#    
#    nx = volume.get_xsize()
#    projections = [] 
#    for iProj in xrange(len(proj_directions)):
#        proj_directions[iProj][2] = random() * 360.0
#        proj_directions[iProj][3] = randrange(-2, 3)
#        proj_directions[iProj][4] = randrange(-2, 3)
#        if random() < 0.5:
#            proj_directions[iProj] = mirror_projection_direction(proj_directions[iProj])
#        projections.append( project(volume, proj_directions[iProj], nx // 2 - 1) )
#        projections[iProj] = filt_tanl(projections[iProj], tanl_filter_cutoff, tanl_filter_fall_off)
#        
#    if output_filename != None and output_filename != "":
#        for iProj in xrange(len(projections)):
#            projections[iProj].write_image(output_filename, iProj)
#        
#    return projections


def get_projections_parameters(projections):
    """
    Read projections parameters from stack of projections
    Returns: list of parameters (phi, theta, psi, sx, sy)
    """
    from EMAN2 import EMData
    import types
    
    if type(projections) is types.StringType:
        projections = EMData.read_images(projections)    
        
    results = []
    for p in projections:
        if not p.has_attr("xform.projection"):
            return None
        param = []
        t = p.get_attr("xform.projection")
        param.append( t.get_params("spider")["phi"  ] )
        param.append( t.get_params("spider")["theta"] )
        param.append( t.get_params("spider")["psi"  ] )
        param.append( t.get_params("spider")["tx"   ] )
        param.append( t.get_params("spider")["ty"   ] )
        results.append(param)
    return results


def angle_between_points_on_sphere(position1, position2):
    """
    Calculate angle>=0 between points on sphere, positions may be lists or transforms
    Returns: angle in degrees
    """
    from math import acos, pi
    from EMAN2 import Transform
    import types
    
    if type(position1) is types.ListType:
        if len(position1) > 2:
            psi = position1[2]
        else:
            psi = 0.0
        position1 = Transform({"type":"spider", "phi":position1[0], "theta":position1[1], "psi":psi})
    
    if type(position2) is types.ListType:
        if len(position2) > 2:
            psi = position2[2]
        else:
            psi = 0.0
        position2 = Transform({"type":"spider", "phi":position2[0], "theta":position2[1], "psi":psi})
   
    unit_vector = [0.0, 0.0, 1.0]
    v1 = position1.inverse() * unit_vector
    v2 = position2.inverse() * unit_vector
    val = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]
    val = max( -1.0, min(1.0, val) )
    return ( acos(val) * (180.0 / pi) )


def choose_nearest_neighbours_from_discrete_angular_space(positions, discrete_angular_space):
    """
    Choose nearest neighbors of given positions from discrete angular space
    Positions is list of lists or transforms
    Returns: list of indexes of points from discrete angular space 
    """
    from EMAN2 import Transform
    import types
    
    result = []
    for p in positions:
        
        if type(p) is types.ListType:
            if len(p) > 2:
                psi = p[2]
            else:
                psi = 0.0
            p = Transform({"type":"spider", "phi":p[0], "theta":p[1], "psi":psi})        
    
        min_angle = 360.0
        best_index = -1
        for i in xrange(len(discrete_angular_space)):
            angle = angle_between_points_on_sphere(p, discrete_angular_space[i])
            if angle < min_angle:
                min_angle = angle
                best_index = i
                
        result.append(best_index)
        
    return result


def fitness_function(projections, discrete_angular_space=None, indexes_of_positions=None):
    """
    Calculate fitness function for projections and their positions in discrete angular space
    Warning: projections given as lists are normalized and their projections parameters are modified
    Returns: average value of L1 distance between input and reconstructed projections, the same for L2 distance
    """
    from EMAN2 import EMData
    from utilities import model_circle, set_params_proj
    from statistics import ave_var
    from reconstruction import recons3d_wbp
    from projection import project
    from global_def import Util
    import types

    if type(projections) is types.StringType:
        projections = EMData.read_images(projections)

    if discrete_angular_space == None:
        discrete_angular_space = get_projections_parameters(projections)
        
    if indexes_of_positions == None:
        indexes_of_positions = range(len(projections))

    nx = projections[0].get_xsize()
    radius = nx // 2 - 1
    mask = model_circle(radius, nx, nx)
    
    avg_proj = ave_var(projections, "")[0]
    avg_value = Util.infomask(avg_proj, mask, True)[0]
    for iProj in xrange(len(projections)):
        p = projections[iProj]
        Util.mul_scalar(p, 1.0 / avg_value)
        set_params_proj(p, discrete_angular_space[indexes_of_positions[iProj]])
    
    vol3d = recons3d_wbp( projections, range(len(projections)) )
    Util.mul_scalar(vol3d, 1.0 / len(projections))
    
    error_L1 = 0.0
    error_L2 = 0.0
    for iProj in xrange(len(projections)):
        p = project(vol3d, discrete_angular_space[indexes_of_positions[iProj]], radius)
        Util.mul_scalar(p, 1.0 / Util.infomask(p, mask, True)[0])
        e1 = projections[iProj].cmp("lod"        ,p,{"mask":mask,"negative":0,"normalize":0}) * (-2.0)
        e2 = projections[iProj].cmp("SqEuclidean",p,{"mask":mask,"zeromask":0,"normto":0})
        error_L1 += e1
        error_L2 += e2
        print "e1/e2=", e1, e2
        
    print "Error (L1/L2)=", error_L1, error_L2
    return error_L1, error_L2, vol3d


def remove_psi(p):
    from EMAN2 import Transform
    from math import pi, asin, acos, cos
    t = Transform({"type":"spider", "phi":p[0], "theta":p[1], "psi":p[2]})
    u = [1.0, 0.0, 0.0]
    v = t * u
    phi = asin( max(-1.0, min(1.0, -v[1])) )
    if v[0] < 0:
        phi = pi - phi
    theta_sin = asin( max(-1.0, min(1.0, v[2] / cos(phi))) )
    theta     = acos( max(-1.0, min(1.0, v[0] / cos(phi))) )
    if (theta_sin < 0):
        theta = 2*pi - theta
    p[0] = phi / pi * 180.0
    p[1] = theta / pi * 180.0
    p[2] = 0.0
    

def rotate_points_on_sphere(points, phi, theta):
    from EMAN2 import Transform
    tr = Transform({"type":"spider", "phi":phi, "theta":theta, "psi":0.0})
    results = []
    for p in points:
        tp = Transform({"type":"spider", "phi":p[0], "theta":p[1], "psi":0.0})
        tp = tr * tp
        p2 = p[:]
        p2[0] = tp.get_params("spider")["phi"]
        p2[1] = tp.get_params("spider")["theta"]
        p2[2] = tp.get_params("spider")["psi"]
        remove_psi(p2)
        results.append(p2)
    return results


def draw_points_on_sphere(points, nx = 64):
    from EMAN2 import Transform
    from utilities import model_blank

    vol = model_blank(nx, nx, nx)
    r = float(nx // 2 - 1)
    u = [0.0, 0.0, 1.0]
    for p in points:
        psi = 0.0
        tr = Transform({"type":"spider", "phi":p[0], "theta":p[1], "psi":psi})
        v = tr.inverse() * u
        x = int(v[0]*r) + nx//2
        y = int(v[1]*r) + nx//2
        z = int(v[2]*r) + nx//2
        vol.set_value_at(x,y,z, 1.0)
        
    return vol


def take_elements(list_with_elements, list_with_indexes):
    result = []
    for i in list_with_indexes:
        result.append(list_with_elements[i])
    return result

 
# returns: error angle for each projection, new transform for each projection 
def quality_of_solution(found_projs_params, original_projs_params, radius=31):
    from utilities import rotation_between_anglesets
    from EMAN2 import Transform
    from pixel_error import max_3D_pixel_error
    import types

    n = len(original_projs_params)

    if not (type(found_projs_params[0]) is types.ListType):
        # convert Transforms to lists
        for i in xrange(n):
            t = [ 1, 2, 3, 0.0, 0.0 ]
            t[0] = found_projs_params[i].get_params("spider")["phi"]
            t[1] = found_projs_params[i].get_params("spider")["theta"]
            t[2] = found_projs_params[i].get_params("spider")["psi"]
            found_projs_params[i] = t

    if type(original_projs_params[0]) is types.ListType:
        # convert lists to Transforms
        for i in xrange(n):
            t = Transform({"type":"spider", "phi":original_projs_params[i][0], "theta":original_projs_params[i][1], "psi":original_projs_params[i][2]})
            original_projs_params[i] = t


    # calculate mirror of given solution 
    found_projs_params_mirror = []
    for p in found_projs_params:
        pm = p[:]
        pm[2] += 180.0
        if (pm[2] >= 360.0):
            pm[2] -= 360.0; 
        found_projs_params_mirror.append(pm)
    
    # rotate solution
    phi, theta, psi = rotation_between_anglesets(found_projs_params, original_projs_params)
#    print "rba=", phi, theta, psi
    tr = Transform({"type":"spider", "phi":phi, "theta":theta, "psi":psi})
    rotated_proj = []
    for p in found_projs_params:
        tp = Transform({"type":"spider", "phi":p[0], "theta":p[1], "psi":p[2]})
        tp = tp * tr
        rotated_proj.append(tp)

    # rotate mirror of solution
    phi, theta, psi = rotation_between_anglesets(found_projs_params_mirror, original_projs_params)
#    print "rba=", phi, theta, psi
    tr = Transform({"type":"spider", "phi":phi, "theta":theta, "psi":psi})
    rotated_proj_mirror = []
    for p in found_projs_params_mirror:
        tp = Transform({"type":"spider", "phi":p[0], "theta":p[1], "psi":p[2]})
        tp = tp * tr
        rotated_proj_mirror.append(tp)
    
    diff = []
    diff_mirror = []
    diff_psi = []
    diff_psi_mirror = []
    
    for i in xrange(n):
        diff           .append( angle_between_points_on_sphere(original_projs_params[i], rotated_proj       [i]) )
        diff_mirror    .append( angle_between_points_on_sphere(original_projs_params[i], rotated_proj_mirror[i]) )
        diff_psi       .append( abs(original_projs_params[i].get_params("spider")["psi"] - rotated_proj       [i].get_params("spider")["psi"]) )
        diff_psi_mirror.append( abs(original_projs_params[i].get_params("spider")["psi"] - rotated_proj_mirror[i].get_params("spider")["psi"]) )
        if diff_psi[i] > 180.0:
            diff_psi[i] = 360.0 - diff_psi[i]
        if diff_psi_mirror[i] > 180.0:
            diff_psi_mirror[i] = 360.0 - diff_psi_mirror[i]        
    
#    print "avg(angle):", sum(diff)/len(diff), sum(diff_mirror)/len(diff_mirror) 
    
    if sum(diff_psi) > sum(diff_psi_mirror):
        diff = diff_mirror
        rotated_proj = rotated_proj_mirror
        diff_psi = diff_psi_mirror
        
    return diff, rotated_proj, diff_psi
    
    # calculate Max Pixel Error for solution and mirror of solution
#    mpe = []
#    mpe_mirror = []
#    for i in xrange(n):
#        mpe       .append( max_3D_pixel_error(original_projs_params[i], rotated_proj       [i], radius) )
#        mpe_mirror.append( max_3D_pixel_error(original_projs_params[i], rotated_proj_mirror[i], radius) )
#        
#    # choose better one from solution and its mirror 
#    print "sum(MPE):", sum(mpe), sum(mpe_mirror)
#    if sum(mpe) > sum(mpe_mirror):
#        mpe = mpe_mirror
#        rotated_proj = rotated_proj_mirror
#        found_projs_params = found_projs_params_mirror
#    del mpe_mirror, rotated_proj_mirror, found_projs_params_mirror
#    
#    # return MPE    
#    return mpe


def mean_stddev_min_max(list_of_elements):
    from math import sqrt
    s = sum(list_of_elements)
    n = len(list_of_elements)
    mean = s / n
    s = 0.0
    min_e = list_of_elements[0]
    max_e = list_of_elements[0]
    for e in list_of_elements:
        s += (e - mean) * (e - mean)
        min_e = min(e, min_e)
        max_e = max(e, max_e)
    if n < 2:
        stddev = 0.0
    else:
        stddev = sqrt(s / (n-1))
    return mean, stddev, min_e, max_e
        

def test_recons3d_wbp(list_proj):
    from EMAN2 import EMData
    from global_def import Util
    from utilities import model_circle, model_blank
    
    B = list_proj[0]
    ny = B.get_ysize()  # have to take ysize, because xsize is different for real and fft images
    nimages = len(list_proj)
    radius = ny // 2 - 1
    mask = model_circle(radius, ny, ny, ny)
    mask2D = model_circle(radius, ny, ny)

    for i in xrange(nimages):
        Util.mul_scalar( list_proj[i], 1.0 / Util.infomask(list_proj[i], mask2D, True)[0] )


    CUBE = EMData()
    CUBE.set_size(ny, ny, ny)
    CUBE.to_zero()

    ss = []
    for iProj in xrange(nimages):
        B = list_proj[iProj]
        transform = B.get_attr("xform.projection")
        d = transform.get_params("spider")
        DMnSS = Util.CANG(d["phi"], d["theta"], d["psi"])
        ss[ (iProj)*6 : (iProj+1)*6 ] = DMnSS["SS"]

    const = 1.0E4

    for iProj in xrange(nimages):
        proj = list_proj[iProj]
        B = proj.copy()
        oldCUBE = CUBE.copy()
        Util.WTF(B, ss, const, iProj+1)  # counting in WTF start from 1!
        Util.BPCQ(B, CUBE, radius)
        diffCUBE = CUBE.copy()
        Util.sub_img(diffCUBE, oldCUBE)
        newCUBE = EMData()
        newCUBE.set_size(ny, ny, ny)
        newCUBE.to_zero()
        Util.BPCQ(B, newCUBE, radius)
        print "error=", newCUBE.cmp("lod"        ,diffCUBE,{"mask":mask,"negative":0,"normalize":0}) * (-2.0)
        print "whole=", newCUBE.cmp("lod"        ,model_blank(ny,ny,ny),{"mask":mask,"negative":0,"normalize":0}) * (-2.0)
        print "errorSqEq=", newCUBE.cmp("SqEuclidean",diffCUBE,{"mask":mask,"zeromask":0,"normto":0})
