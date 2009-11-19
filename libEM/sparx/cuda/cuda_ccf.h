void calculate_ccf(float *subject_image, float *ref_image, float *ccf, int NIAMGE, int NX, int NY, int RING_LENGTH, int NRING, int OU, float STEP, int KX, int KY, float *sx, float *sy, int id, int silent); 

void filter_image(float *image_in, float *image_out, int NIMA, int NX, int NY, float *params, int id);

void rot_filt_sum(float *image, int NIMA, int NX, int NY, int CTF, float *ctf_params, float *ali_params, float *ave1, float *ave2, int id);
