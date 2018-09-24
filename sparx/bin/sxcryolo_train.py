#
# Copyright (C) 2016  Thorsten Wagner (thorsten.wagner@mpi-dortmund.mpg.de)
#

# source code in this file under either license. However, note that the
# complete EMAN2 and SPARX software packages have some GPL dependencies,
# so you are responsible for compliance with the licenses of these packages
# if you opt to use BSD licensing. The warranty disclaimer below holds
# in either instance.
#
# This complete copyright notice must be included in any revised version of the
# source code. Additional authorship citations may be added, but existing
# author citations must be preserved.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

from __future__ import print_function
import argparse
from json import dump
import subprocess

argparser = argparse.ArgumentParser(
    description='Train crYOLO on any dataset',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

argparser.add_argument(
    'architecture',
    type=str,
    help='Specifiy the used model')

argparser.add_argument(
    'particle_diameter',
    type=int,
    help='Particle diameter in pixel')

argparser.add_argument(
    'training_dir',
    type=str,
    help='Path to training images')

argparser.add_argument(
    'annot_dir',
    type=str,
    help='Path to training images')

argparser.add_argument(
    '--input_size',
    type=int,
    default=768,
    help='Size to which the image is resized')

argparser.add_argument(
    '--num_patches',
    type=int,
    default=2,
    help='Number of patches used during training')

argparser.add_argument(
    '--train_times',
    type=int,
    default=10,
    help='How often a images is augmented and repeadet in one epoch.')

argparser.add_argument(
    '--weights_name',
    type=str,
    default="cryolo_model.h5",
    help='Name of the model')

argparser.add_argument(
    '--batch_size',
    type=int,
    default=6,
    help='How many patches are processed in parallel')


argparser.add_argument(
    '--learning_rate',
    type=float,
    default=0.0001,
    help='Learning rate used during training.')

argparser.add_argument(
    '--np_epoch',
    type=int,
    default=100,
    help='Maximum number of epochs.')

argparser.add_argument(
    '--object_scale',
    type=float,
    default=5.0,
    help='Loss scale for particles.')

argparser.add_argument(
    '--no_object_scale',
    type=float,
    default=1.0,
    help='Loss scale for background.')

argparser.add_argument(
    '--coord_scale',
    type=float,
    default=1.0,
    help='Loss scale for coordinates.')

argparser.add_argument(
    '--valid_image_dir',
    type=str,
    help='Path to validation images')

argparser.add_argument(
    '--valid_annot_dir',
    type=str,
    help='Path to validation annotations')

argparser.add_argument(
    '--warmup',
    type=int,
    help='Warm up epochs')

argparser.add_argument(
    '--gpu',
    default=0,
    type=int,
    nargs="+",
    help="Specifiy which gpu(s) should be used. Multiple GPUs are separated by a whitespace")

argparser.add_argument(
    '--early',
    default=5,
    type=int,
    help='Early stop patience. If the validation loss did not improve longer than the early stop patience, '
         'the training is stopped.')


def main():
    #Read arguments
    args = argparser.parse_args()

    architecture = args.architecture
    particle_diameter = args.particle_diameter
    trainging_dir = args.training_dir
    annot_dir = args.annot_dir
    input_size = args.input_size
    num_patches = args.num_patches
    train_times = args.train_times
    weights_name = args.weights_name
    batch_size = args.batch_size
    learning_rate = args.learning_rate
    np_epoch = args.np_epoch
    object_scale = args.object_scale
    no_object_scale = args.no_object_scale
    coord_scale = args.coord_scale
    valid_image_dir = args.valid_image_dir
    valid_annot_dir = args.valid_annot_dir
    warmup = args.warmup
    early_stop = int(args.early)

    if type(args.gpu) is list:
        str_gpus = [str(entry) for entry in args.gpu]
    else:
        str_gpus = str(args.gpu)

    arg_gpu = ' '.join(str_gpus)

    #Create config file
    model_dict = {'architecture': architecture,
                  'input_size': input_size,
                  'anchors': [particle_diameter, particle_diameter],
                  'overlap_patches': 200,
                  'max_box_per_image': 700,
                  'num_patches': num_patches}

    train_dict = {'train_image_folder': trainging_dir,
                  'train_annot_folder': annot_dir,
                  'train_times': train_times,
                  'pretrained_weights': weights_name,
                  'batch_size': batch_size,
                  'learning_rate': learning_rate,
                  'nb_epoch': np_epoch,
                  'warmup_epochs': 0,
                  'object_scale': object_scale,
                  'no_object_scale': no_object_scale,
                  'coord_scale': coord_scale,
                  'class_scale': 1.0,
                  "saved_weights_name": weights_name,
                  "debug": True,
                  "log_path": "cryolo_logs/"
                  }

    valid_dict = {'valid_image_folder': valid_image_dir,
                  'valid_annot_folder': valid_annot_dir,
                  'valid_times': 1
                  }
    dict = {"model": model_dict, "train": train_dict, "valid": valid_dict}
    path = "config_yolo.json"
    with open(path, 'w') as f:
        dump(dict, f, ensure_ascii=False, indent=4)

    #Run the training
    warmup_argument = "-w=" + str(warmup)
    gpu_argument = "-g=" + arg_gpu
    early_stop = "-e=" + str(early_stop)
    subprocess.check_call(['python', 'cryolo_train.py', "-c=config_yolo.json", warmup_argument, gpu_argument, early_stop])
    warmup_argument = "-w=" + 0
    subprocess.check_call(['python', 'cryolo_train.py', "-c=config_yolo.json", warmup_argument, gpu_argument, early_stop])





