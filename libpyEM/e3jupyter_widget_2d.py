#!/usr/bin/env python
#
# Author: Anya Porter, 06/11/2024 (anastasia.porter@bcm.edu)
# Copyright (c) 2000-2006 Baylor College of Medicine
#
# This software is issued under a joint BSD/GNU license. You may use the
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
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston MA 02111-1307 USA
#
#

from EMAN2 import *
import numpy as np
import matplotlib.pyplot as plt
import ipywidgets as widget
from ipyevents import Event
from IPython.display import display
import asyncio
from time import time
from io import BytesIO
from PIL import Image, ImageDraw

def controls_layout(above):
	"""Sets the layout for the controls box"""
	if above:
		return widget.Layout(width='100%', border='solid 1px black', margin='0px 10px 10px 0px', padding='5px 5px 5px 5px')
	else:
		return widget.Layout(width='36%', border='solid 1px black', margin='0px 10px 10px 0px', padding='5px 5px 5px 5px')

def image_layout(above):
	"""Sets the layout for the image box"""
	if above:
		return widget.Layout(border='solid 1px black', margin='0px 10px 10px 0px', object_fit='none', object_position='0px 0px')
	else:
		return widget.Layout(width='64%', border='solid 1px black', margin='0px 10px 10px 0px', object_fit='none', object_position='0px 0px')

def button_layout(above):
	"""Sets the layout for buttons"""
	if above:
		return widget.Layout(width='100%')
	else:
		return widget.Layout(width='32%')

def slider_layout():
	"""Sets the layout for sliders"""
	return widget.Layout(width="99%")

def text_layout():
	"""Sets the layout for text"""
	return widget.Layout(justify_content='space-around', align_content='flex-start')


def normalize(numpy_eman, min, max, invert=False):
	"""Takes in numpy_eman (a numpy float array), a minimum value, a maximimum value, and if it should be inverted
	and returns an 8bit numpy array that is that mapping"""
	if min == max:
		max += 0.001
	if invert:
		return ((np.clip(numpy_eman, min, max) - max)*255/(min-max)).astype(np.uint8)
	else:
		return ((np.clip(numpy_eman, min, max) - min)*255/(max - min)).astype(np.uint8)

def convert_to_bytes(img_np, scale, mode=0, xy=[], color=False):
	"""Takes in an 8bit numpy array, a scaling factor, a mode (0=app, 1=probe, 2=measure), a set of coordinates xy,
	and a bool for the color and returns a byte representation of the image at the given scale with any probe/measure lines drawn"""
	in_mem_file = BytesIO()
	img_pil = Image.fromarray(img_np, 'L')
	img_pil = img_pil.transpose(Image.FLIP_TOP_BOTTOM)
	if mode == 1 and len(xy) == 2: # Probe
		if color:
			outline_color = 0
		else:
			outline_color = 255
		draw = ImageDraw.Draw(img_pil)
		draw.rectangle(xy, fill=None, outline=outline_color, width=max(1, int(1/scale)))
	if mode == 2 and len(xy) == 2: # Meas
		if color:
			outline_color = 0
		else:
			outline_color = 255
		draw = ImageDraw.Draw(img_pil)
		draw.line(xy, fill=outline_color, width=max(1, int(1/scale)))
	img_pil = img_pil.resize((int(img_np.shape[1]*scale), int(img_np.shape[0]*scale)), resample=Image.Resampling.NEAREST)
	img_pil.save(in_mem_file, format="PNG")
	in_mem_file.seek(0)
	return in_mem_file.read()

################### Code from Ipywidgets documentation ######################
class Timer:
	def __init__(self, timeout, callback):
		self._timeout = timeout
		self._callback = callback

	async def _job(self):
		await asyncio.sleep(self._timeout)
		self._callback()

	def start(self):
		self._task = asyncio.ensure_future(self._job())

	def cancel(self):
		self._task.cancel()

def debounce(wait):
	""" Decorator that will postpone a function's
		execution until after `wait` seconds
		have elapsed since the last time it was invoked. """
	def decorator(fn):
		timer = None
		def debounced(*args, **kwargs):
			nonlocal timer
			def call_it():
				fn(*args, **kwargs)
			if timer is not None:
				timer.cancel()
			timer = Timer(wait, call_it)
			timer.start()
		return debounced
	return decorator

def throttle(wait):
	""" Decorator that prevents a function from being called
		more than once every wait period. """
	def decorator(fn):
		time_of_last_call = 0
		scheduled, timer = False, None
		new_args, new_kwargs = None, None
		def throttled(*args, **kwargs):
			nonlocal new_args, new_kwargs, time_of_last_call, scheduled, timer
			def call_it():
				nonlocal new_args, new_kwargs, time_of_last_call, scheduled, timer
				time_of_last_call = time()
				fn(*new_args, **new_kwargs)
				scheduled = False
			time_since_last_call = time() - time_of_last_call
			new_args, new_kwargs = args, kwargs
			if not scheduled:
				scheduled = True
				new_wait = max(0, wait - time_since_last_call)
				timer = Timer(new_wait, call_it)
				timer.start()
		return throttled
	return decorator

################### End code from Ipywidgets documentation ######################
class JupyterDisplay():
	"""Designed to display an Eman object in Jupyter with a control panel"""
	output = widget.Output()

	def __init__(self, eman_img_stack, controls_above=False):
		"""Initialize the object and create the plot and widgets
		Inputs:
			self
			eman_img_stack - a list of EmanData objects, one that "Show2d" would be an option in e2display.py
		Outputs:
			none, displays widget when called"""

		style={'description_width': 'initial', 'handle_color': 'grey'}

		# Define widgets
		self.controls_above = controls_above
		if isinstance(eman_img_stack, list):
			self.eman_object = eman_img_stack
		else:
			self.eman_object = [eman_img_stack]
		self.index = widget.IntSlider(min=0, max=len(self.eman_object)-1,value=0, description='N#',
						layout=slider_layout(), style=style)
		self.cur_image = self.eman_object[self.index.value].numpy()
		self.fft = None
		self.fft_display = None
		self.image_min = self.eman_object[self.index.value]['minimum']
		self.image_max = self.eman_object[self.index.value]['maximum']
		self.mean = self.eman_object[self.index.value]['mean']
		self.sigma = self.eman_object[self.index.value]['sigma']
		self.fourier_min = None
		self.fourier_max = None
		self.fourier_mean = None
		self.fourier_sigma = None
		self.min_slider = widget.FloatSlider(min=self.image_min, max=self.image_max, step=0.001, value=self.image_min, description="Min",
							readout_format=".3f", layout=slider_layout(), style=style)
		self.max_slider = widget.FloatSlider(min=self.image_min, max=self.image_max, step=0.001, value=self.image_max, description="Max",
							readout_format=".3f", layout=slider_layout(), style=style)
		self.brightness = widget.FloatSlider(min=-1, max=1, step=0.001, value=0, description="Brt", readout_format=".3f",
							layout=slider_layout(), style=style)
		self.contrast = widget.FloatSlider(min=0, max=1, step=0.001, value=0.5, description="Cont", readout_format=".3f",
							layout=slider_layout(), style=style)
		self.slider_tabs = widget.Tab()
		self.slider_tabs.children = [widget.VBox([self.min_slider, self.max_slider]),
						widget.VBox([self.brightness, self.contrast])]
		self.slider_tabs.titles = ["MinMax", "BrtCont"]
		self.auto_button = widget.Button(description="autoC", layout=button_layout(self.controls_above))
		self.full_button = widget.Button(description="fullC", layout=button_layout(self.controls_above))
		self.invert_button = widget.ToggleButton(description="Invert", layout=button_layout(self.controls_above))
		self.lock = False
		self.fourier_buttons = widget.ToggleButtons(options=['Real', 'Amp', 'Pha'])
		if self.controls_above:
			self.fourier_buttons.style={'button_width': '100%'}
		else:
			self.fourier_buttons.style={'button_width': '32%'}
		self.mag_slider = widget.FloatSlider(min=0.1, max=10.0, value=1, description="Mag",
							layout=slider_layout(), style=style)
		# App tab
		x,y,z = self.eman_object[self.index.value].get_sizes()
		self.apix_x = self.eman_object[self.index.value]['apix_x']
		self.apix_y = self.eman_object[self.index.value]['apix_y']
		self.apix_z = self.eman_object[self.index.value]['apix_z']
		try:
			ctf_defocus = self.eman_object[self.index.value]['ctf'].defocus
			delta_unicode = "\u0394"
			size = widget.VBox([widget.Label(value=f'nx={x:d}'),
							widget.Label(value=f'ny={y:d}'),
							widget.Label(value=f'nz={z:d}'),
							widget.Label(value=f'{delta_unicode}Z={ctf_defocus:1.5g}')],
							layout=text_layout())
		except:
			size = widget.VBox([widget.Label(value=f'nx={x:d}'),
							widget.Label(value=f'nx={y:d}'),
							widget.Label(value=f'nx={z:d}')],
							layout=text_layout())
		stats = widget.VBox([widget.Label(value=f'min={self.image_min:1.4g}'),
							widget.Label(value=f'max={self.image_max:1.4g}'),
							widget.Label(value=f'mean={self.mean:1.4g}'),
							widget.Label(value=f'sigma={self.sigma:1.5g}')],
							layout=text_layout())
		try:
			particle_repr = self.eman_object[self.index.value]['ptcl_repr']
			apix = widget.VBox([widget.Label(value=f'apix_x={self.apix_x:1.3f}'),
							widget.Label(value=f'apix_y={self.apix_y:1.3f}'),
							widget.Label(value=f'apix_z={self.apix_z:1.3f}'),
							widget.Label(value=f'ptcl_repr={particle_repr}')],
							layout=text_layout())
		except:
			apix = widget.VBox([widget.Label(value=f'apix_x={self.apix_x:1.3f}'),
							widget.Label(value=f'apix_y={self.apix_y:1.3f}'),
							widget.Label(value=f'apix_z={self.apix_z:1.3f}')],
							layout=text_layout())
		app_tab = widget.HBox([size, stats, apix], layout=text_layout())
		# Probe tab
		self.line_color_button = widget.ToggleButton(description="Line Col")
		self.probe_size_widget = widget.IntText(value=32, description="Probe Size")
		left_probe_col = widget.VBox([widget.Label(value="Area Avg: "),
							widget.Label(value="Area Sig: "),
							widget.Label(value="Skewness: "),
							widget.Label(value="Kurtosis: ")],
							layout=text_layout())
		right_probe_col = widget.VBox([widget.Label(value="Area Avg (!=0): "),
							widget.Label(value="Area Sig (!=0): "),
							widget.Label(value="Center Coord: "),
							widget.Label(value="dcen () ")],
							layout=text_layout())
		probe_tab = widget.VBox([self.probe_size_widget,
							widget.HBox([widget.Label(value= "Point Value (ctr pix): "), self.line_color_button]),
							widget.HBox([left_probe_col, right_probe_col])], layout=text_layout())
		# Measure tab
		self.aPix_slider = widget.FloatSlider(value=self.apix_x , min=0.5, max=10, step=0.001, description="A/Pix", 
							readout_format=".3f", layout=slider_layout(), style=style)
		left_meas_col = widget.VBox([widget.Label(value="Start: 0,0"),
							widget.Label(value="dx,dx: 0"),
							widget.Label(value="Value: ?")], layout=text_layout())
		right_meas_col = widget.VBox([widget.Label(value="End 0,0"),
							widget.Label(value="Len: 0")], layout=text_layout())
		meas_tab = widget.VBox([widget.HBox([self.aPix_slider, self.line_color_button]),
							widget.HBox([left_meas_col, right_meas_col], layout=text_layout())])
		self.tabs = widget.Tab()
		self.tabs.children = [app_tab, probe_tab, meas_tab]
		self.tabs.titles = ["App", "Probe", "Meas"]
		if self.controls_above:
			self.tabs.layout = widget.Layout(width='35%')
		self.start_coord = None
		self.coords = []

		self.img_norm = normalize(self.cur_image, self.min_slider.value, self.max_slider.value, self.invert_button.value)
		self.img_widget = widget.Image(value=convert_to_bytes(self.img_norm, self.mag_slider.value), format='png')


		# Create the histogram
		with plt.ioff():
			self.hist_fig, self.hist_axs = plt.subplots(1, 1, figsize=(3, 2))
			self.hist_axs.hist(self.cur_image.flatten(), bins='auto')
			self.hist_fig.canvas.toolbar_position = "bottom"
		plt.tight_layout()

		# Observe changes
		self.invert_button.observe(self.invert, 'value')
		self.auto_button.on_click(self.autoC)
		self.full_button.on_click(self.fullC)
		self.brightness.observe(self.update_brt, 'value')
		self.contrast.observe(self.update_cont, 'value')
		self.min_slider.observe(self.update_min, 'value')
		self.max_slider.observe(self.update_max, 'value')
		self.index.observe(self.change_image, 'value')
		self.index.observe(self.update_app_stats, 'value')
		self.fourier_buttons.observe(self.fourier_transform, 'value')
		# self.tabs.observe(self.set_mode, 'selected_index')
		self.mag_slider.observe(self.update_zoom_slider, 'value')

		# Start with auto contrast
		self.auto_button.click()

		# Define all events (won't be connected unless called)
		self.no_drag = Event(source=self.img_widget, watched_events=['dragstart'], prevent_default_action=True)
		self.mouse_updown = Event(source=self.img_widget, watched_events=['mousedown', 'mouseup'])
		self.mouse_move = Event(source=self.img_widget, watched_events=['mousemove'], wait=500)
		self.mouse_scroll = Event(source=self.img_widget, watched_events=['wheel'], prevent_default_action=True)
		self.control_position = Event(source=self.img_widget, watched_events=['auxclick'], prevent_default_action=True)
		self.mouse_updown.on_dom_event(self.connect_move_event)
		self.mouse_scroll.on_dom_event(self.update_zoom_scroll)
		self.control_position.on_dom_event(self.display_control_panel)



		# Establish Layout
		with self.output:
			if self.controls_above:
				contrast_box = widget.VBox([self.fourier_buttons, self.invert_button, self.auto_button, self.full_button],
								layout=widget.Layout(width='10%'))
				self.tabs.layout = widget.Layout(width='35%')
				self.fourier_buttons.style={'button_width': '99%'}
				# sliders = widget.VBox([self.mag_slider, self.min_slider, self.max_slider, self.brightness,
				#						self.contrast, self.index], layout=widget.Layout(width='33%'))
				sliders = widget.VBox([self.mag_slider, self.slider_tabs, self.index], layout=widget.Layout(width='33%'))
				controls = widget.HBox([self.tabs, contrast_box, sliders, self.hist_fig.canvas])
				controls.layout = controls_layout(self.controls_above)
				self.img_widget.layout = image_layout(self.controls_above)
				display(widget.VBox([controls, self.img_widget]))
			else:
				contrast_box = widget.HBox([self.invert_button, self.auto_button, self.full_button])
				# controls = widget.VBox([self.tabs, self.fourier_buttons, contrast_box, self.mag_slider,
				#						 self.min_slider, self.max_slider, self.brightness, self.contrast,
				#						 self.hist_fig.canvas, self.index])
				self.tabs.layout = widget.Layout(width='100%')
				self.slider_tabs.layout = widget.Layout(width='100%')
				self.fourier_buttons.style={'button_width': '32%'}
				controls = widget.VBox([self.tabs, self.fourier_buttons, contrast_box, self.mag_slider,
								self.slider_tabs, self.hist_fig.canvas, self.index])
				controls.layout = controls_layout(self.controls_above)
				self.img_widget.layout = image_layout(self.controls_above)
				display(widget.HBox([controls, self.img_widget]))
		display(self.output)


	def __del__(self):
		"""When object is deleted tell jupyter to close the figures"""
		plt.close(self.hist)




	@output.capture()
	def draw_image(self):
		if self.slider_tabs.selected_index == 0:
			# Min/max mode
			if self.fourier_buttons.value == 'Real':
				self.img_norm = normalize(self.cur_image, self.min_slider.value, self.max_slider.value, self.invert_button.value)
			else:
				self.img_norm = normalize(self.fft_display, self.min_slider.value, self.max_slider.value, self.invert_button.value)
		elif self.slider_tabs.selected_index == 1:
			# Brt/cont mode
			if self.fourier_buttons.value == 'Real':
				min = ((self.image_min+self.image_max) / 2.0)-(self.image_max-self.image_min)*(1.0-self.contrast.value)-self.brightness.value*(self.image_max-self.image_min)
				max = ((self.image_min+self.image_max) / 2.0)+(self.image_max-self.image_min)*(1.0-self.contrast.value)-self.brightness.value*(self.image_max-self.image_min)
				self.img_norm = normalize(self.cur_image, min, max, self.invert_button.value)
			else:
				min = ((self.fourier_min+self.fourier_max) / 2.0)-(self.fourier_max-self.fourier_min)*(1.0-cont)-brt*(self.fourier_max-self.fourier_min)
				max = ((self.fourier_min+self.fourier_max) / 2.0)+(self.fourier_max-self.fourier_min)*(1.0-cont)-brt*(self.fourier_max-self.fourier_min)
				self.img_norm = normalize(self.fft_display, min, max, self.invert_button.value)
		self.img_widget.value = convert_to_bytes(self.img_norm, self.mag_slider.value, mode=self.tabs.selected_index, xy=self.coords, color=self.line_color_button.value)

	def update_slider_mode(self, change):
		# 0 = min/max, 1 = brt/cont
		if change['new'] == 0:
			self.update_min_max(self.brightness.value, self.contrast.value)
			self.brightness.unobserve(self.update_brt, 'value')
			self.contrast.unobserve(self.update_cont, 'value')
			self.min_slider.observe(self.update_min, 'value')
			self.max_slider.observe(self.update_max, 'value')
		elif change['new'] == 1:
			self.update_brt_cont(self.min_slider.value, self.max_slider.value)
			self.brightness.observe(self.update_brt, 'value')
			self.contrast.observe(self.update_cont, 'value')
			self.min_slider.unobserve(self.update_min, 'value')
			self.max_slider.unobserve(self.update_max, 'value')

	# Functions for widgets
	def update_app_stats(self, change):
		"""Updates the App tab based off a change in self.index changing which
		image in the stack is selected."""
		selected_tab = self.tabs.selected_index
		x,y,z = self.eman_object[change['new']].get_sizes()
		self.apix_x = self.eman_object[change['new']]['apix_x']
		self.apix_y = self.eman_object[change['new']]['apix_y']
		self.apix_z = self.eman_object[change['new']]['apix_z']
		try:
			ctf_defocus = self.eman_object[change['new']]['ctf'].defocus
			delta_unicode = "\u0394"
			size = widget.VBox([widget.Label(value=f'nx={x:d}'),
							widget.Label(value=f'ny={y:d}'),
							widget.Label(value=f'nz={z:d}'),
							widget.Label(value=f'{delta_unicode}Z={ctf_defocus:1.5g}')],
							layout=text_layout())
		except:
			size = widget.VBox([widget.Label(value=f'nx={x:d}'),
							widget.Label(value=f'nx={y:d}'),
							widget.Label(value=f'nx={z:d}')],
							layout=text_layout())
		stats = widget.VBox([widget.Label(value=f'min={self.image_min:1.4g}'),
							widget.Label(value=f'max={self.image_max:1.4g}'),
							widget.Label(value=f'mean={self.mean:1.4g}'),
							widget.Label(value=f'sigma={self.sigma:1.5g}')],
							layout=text_layout())
		try:
			particle_repr = self.eman_object[change['new']]['ptcl_repr']
			apix = widget.VBox([widget.Label(value=f'apix_x={self.apix_x:1.3f}'),
							widget.Label(value=f'apix_y={self.apix_y:1.3f}'),
							widget.Label(value=f'apix_z={self.apix_z:1.3f}'),
							widget.Label(value=f'ptcl_repr={particle_repr}')],
							layout=text_layout())
		except:
			apix = widget.VBox([widget.Label(value=f'apix_x={self.apix_x:1.3f}'),
							widget.Label(value=f'apix_y={self.apix_y:1.3f}'),
							widget.Label(value=f'apix_z={self.apix_z:1.3f}')],
							layout=text_layout())
		app_tab = widget.HBox([size, stats, apix], layout=text_layout())
		# with self.tabs.hold_trait_notifications(), self.tabs.hold_sync():
		# with self.tabs._lock_property({'selected_index': selected_tab}):
		if self.tabs.selected_index > 0:
			self.tabs.selected_index = None
		self.tabs.children = [app_tab, self.tabs.children[1], self.tabs.children[2]]
		self.tabs.selected_index = selected_tab

	def change_image(self, change):
		"""Updates the image displayed from the stack of images self.eman_object based on
		the slider self.index. Will keep the previous fourier status. Also updates histogram"""
		self.lock = True
		self.cur_image = self.eman_object[change['new']].numpy()
		self.fft = None
		self.fft_display = None
		if self.fourier_buttons.value == 'Real':
			self.image_min = self.eman_object[change['new']]['minimum']
			self.image_max = self.eman_object[change['new']]['maximum']
			self.mean = self.eman_object[change['new']]['mean']
			self.sigma = self.eman_object[change['new']]['sigma']
			self.hist_axs.cla()
			self.hist_axs.hist(self.cur_image.flatten(), bins='auto')
		else:
			self.fft = self.eman_object[change['new']].do_fft()
			self.fft_display = self.fft.process("xform.fourierorigin.tocenter")
			if self.fourier_buttons.value == 'Amp':
				self.fft_display = self.fft_display.get_fft_amplitude()
			elif self.fourier_buttons.value == 'Pha':
				self.fft_display = self.fft_display.get_fft_phase()
			self.fourier_min = self.fft_display['minimum']
			self.fourier_max = self.fft_display['maximum']
			self.fourier_mean = self.fft_display['mean']
			self.fourier_sigma = self.fft_display['sigma']
			self.fft_display = self.fft_display.numpy().copy()
			self.hist_axs.cla()
			self.hist_axs.hist(self.fft_display.flatten(), bins='auto', range=(self.fourier_min, self.fourier_max))
		self.draw_image()
		self.lock = False

	@output.capture()
	def fourier_transform(self, change):
		"""Switches between displaying the real image, the fourier amplitude, or the fourier phase of the image
		based on the toggle buttons self.fourier_buttons. Also changes histogram."""
		self.lock = True
		if change['new'] == 'Real':
			self.image_min = self.eman_object[self.index.value]['minimum']
			self.image_max = self.eman_object[self.index.value]['maximum']
			self.mean = self.eman_object[self.index.value]['mean']
			self.sigma = self.eman_object[self.index.value]['sigma']
			self.hist_axs.cla()
			self.hist_axs.hist(self.cur_image.flatten(), bins='auto')
		else:
			if self.fft == None:
				self.fft = self.eman_object[self.index.value].do_fft()
			self.fft_display = self.fft.process("xform.fourierorigin.tocenter")
			if change['new'] == 'Amp':
				self.fft_display = self.fft_display.get_fft_amplitude()
			elif change['new'] == 'Pha':
				self.fft_display = self.fft_display.get_fft_phase()
			self.fourier_min = self.fft_display['minimum']
			self.fourier_max = self.fft_display['maximum']
			self.fourier_mean = self.fft_display['mean']
			self.fourier_sigma = self.fft_display['sigma']
			self.fft_display = self.fft_display.numpy().copy()
			self.hist_axs.cla()
			self.hist_axs.hist(self.fft_display.flatten(), bins='auto', range=(self.fourier_min, self.fourier_max))
		self.auto_button.click()
		self.draw_image()
		self.lock = False

	def autoC(self, b):
		"""Sets the self.min_slider and self.max_slider values for auto contrast"""
		status = self.lock
		self.lock = True
		if self.fourier_buttons.value == 'Real':
			self.min_slider.min = self.image_min
			self.min_slider.max = self.image_max
			self.max_slider.min = self.image_min
			self.max_slider.max = self.image_max
			self.min_slider.value = max(self.image_min, self.mean-3.0*self.sigma)
			self.max_slider.value = min(self.image_max, self.mean+3.0*self.sigma)
		else:
			self.min_slider.min = 0
			self.min_slider.max = min(self.fourier_max, self.fourier_mean+20.0*self.fourier_sigma)
			self.max_slider.min = 0
			self.max_slider.max = min(self.fourier_max, self.fourier_mean+20.0*self.fourier_sigma)
			self.min_slider.value = 0
			self.max_slider.value = min(self.fourier_max, self.fourier_mean+4.0*self.fourier_sigma)
		self.update_brt_cont(self.min_slider.value, self.max_slider.value)
		self.lock = status
		if self.lock:
			return
		else:
			if self.fourier_buttons.value == 'Real':
				self.img_norm = normalize(self.cur_image, self.min_slider.value, self.max_slider.value, self.invert_button.value)
			else:
				self.img_norm = normalize(self.fft_display, self.min_slider.value, self.max_slider.value, self.invert_button.value)
			self.img_widget.value = convert_to_bytes(self.img_norm, self.mag_slider.value, mode=self.tabs.selected_index, xy=self.coords, color=self.line_color_button.value)

	def fullC(self, b):
		"""Sets the self.min_slider and self.max_slider values for auto contrast"""
		status = self.lock
		self.lock = True
		if self.fourier_buttons.value == 'Real':
			self.min_slider.min = self.image_min
			self.min_slider.max = self.image_max
			self.max_slider.min = self.image_min
			self.max_slider.max = self.image_max
			self.min_slider.value = self.image_min
			self.max_slider.value = self.image_max
		else:
			self.min_slider.min = 0
			self.max_slider.min = 0
			self.min_slider.max = min(self.fourier_max, self.fourier_mean+20.0*self.fourier_sigma)
			self.max_slider.max = min(self.fourier_max, self.fourier_mean+20.0*self.fourier_sigma)
			self.min_slider.value = 0
			self.max_slider.value = self.fourier_max
		self.update_brt_cont(self.min_slider.value, self.max_slider.value)
		self.lock = status
		if self.lock:
			return
		else:
			self.draw_image()

	def invert(self, change):
		"""Inverts the image"""
		self.draw_image()

	def update_brt_cont(self, min, max):
		"""Updates self.brightness and self.contrast sliders to reflect changes in self.min_slider or
		self.max_slider. Throttled to prevent jumpiness."""
		if self.fourier_buttons.value == 'Real':
			self.brightness.value = -0.5*(min+max-(self.image_max+self.image_min))/((self.image_max-self.image_min))
			self.contrast.value = 1.0 - ((min-max) / (2.0*(self.image_min-self.image_max)))
		else:
			self.brightness.value = -0.5*(min+max-(self.fourier_max+self.fourier_min))/((self.fourier_max-self.fourier_min))
			self.contrast.value = 1.0 - ((min-max) / (2.0*(self.fourier_min-self.fourier_max)))

	def update_min_max(self, brt, cont):
		"""Updates self.min_slider and self.max_slider to reflect changes in self.brightness or
		self.contrast sliders. Throttled to prevent jumpiness."""
		if self.fourier_buttons.value == 'Real':
			m0 = ((self.image_min+self.image_max) / 2.0)-(self.image_max-self.image_min)*(1.0-cont)-brt*(self.image_max-self.image_min)
			m1 = ((self.image_min+self.image_max) / 2.0)+(self.image_max-self.image_min)*(1.0-cont)-brt*(self.image_max-self.image_min)
		else:
			m0 = ((self.fourier_min+self.fourier_max) / 2.0)-(self.fourier_max-self.fourier_min)*(1.0-cont)-brt*(self.fourier_max-self.fourier_min)
			m1 = ((self.fourier_min+self.fourier_max) / 2.0)+(self.fourier_max-self.fourier_min)*(1.0-cont)-brt*(self.fourier_max-self.fourier_min)
		with self.min_slider.hold_trait_notifications():
			self.min_slider.value = m0
			self.min_slider.min = min(m0-self.min_slider.step, self.min_slider.min)
		with self.max_slider.hold_trait_notifications():
			self.max_slider.value = m1
			self.max_slider.max = max(m1+self.max_slider.step, self.max_slider.max)

	def update_min(self, change):
		"""Updates the minimum value mapped to the colormap based on a change to self.min_slider.
		If the change was due to changing min_slider itself it updates brightness and contrast and
		redraws the image."""
		if not self.lock:
			self.lock = True
			# self.update_brt_cont(change['new'], self.max_slider.value)
			self.draw_image()
			self.lock = False

	def update_max(self, change):
		"""Updates the maximum value mapped to the colormap based on a change to self.max_slider.
		If the change was due to changing max_slider itself itself it updates brightness and contrast and
		redraws the image."""
		if not self.lock:
			self.lock = True
			# self.update_brt_cont(self.min_slider.value, change['new'])
			self.draw_image()
			self.lock = False

	def update_brt(self, change):
		"""Updates self.min_slider and self.max_slider when the user changes self.brightness."""
		if not self.lock:
			self.lock = True
			# self.update_min_max(change['new'], self.contrast.value)
			self.draw_image()
			self.lock = False

	def update_cont(self, change):
		"""Updates self.min_slider and self.max_slider when the user changes self.contrast."""
		if not self.lock:
			self.lock = True
			# self.update_min_max(self.brightness.value, change['new'])
			self.draw_image()
			self.lock = False

	def update_zoom_slider(self, change):
		"""Updates the image based on a change in self.mag_slider. A 2x increase in the value will be a 2x increase
		in the image side length"""
		self.img_widget.value = convert_to_bytes(self.img_norm, change['new'], mode=self.tabs.selected_index, xy=self.coords, color=self.line_color_button.value)

	def update_zoom_scroll(self, event):
		"""Updates mag slider (and thus the image) based on a scroll event"""
		# Need to use event['deltaY']--it is negative when scrolling up and positive when scrolling down
		# It is also measured in pixels currently and one move of my mouse is 108 so gives a change of ~0.1
		new_mag = self.mag_slider.value - event['deltaY']/1000
		self.mag_slider.min = max(min(self.mag_slider.min, new_mag), 0.0001)
		self.mag_slider.max = max(self.mag_slider.max, new_mag)
		self.mag_slider.value = new_mag
		# Updates image by causing a change in mag_slider.value

	def pan(self, event):
		"""Moves the image within the box with the mouse. Parts outside the box will get cut off"""
		self.img_widget.layout.object_position = f"{event['relativeX']-self.start_coord[0]}px {event['relativeY']-self.start_coord[1]}px"

	def connect_move_event(self,event):
		"""Based on which tab is selected (App, Probe, or Meas), connects or disconnects
		the mouse_move Event to the proper function"""
		if event['event']=='mousedown':
			# Connect proper events on pressing a mouse button
			if self.tabs.selected_index == 0: # App tab
				original_offset = [int(x[:-2]) for x in self.img_widget.layout.object_position.split()]
				self.start_coord = (event['relativeX'] - original_offset[0], event['relativeY']-original_offset[1])
				self.mouse_move.on_dom_event(self.pan)
			elif self.tabs.selected_index == 1: # Probe tab
				original_offset = [int(x[:-2]) for x in self.img_widget.layout.object_position.split()]
				center = (int((event['relativeX'] - original_offset[0])/self.mag_slider.value),
						int((event['relativeY'] - original_offset[1])/self.mag_slider.value))
				self.coords = [(center[0] - (self.probe_size_widget.value // 2), center[1] - (self.probe_size_widget.value // 2)),
						(center[0] + (self.probe_size_widget.value // 2), center[1] + (self.probe_size_widget.value // 2))]
				self.mouse_move.on_dom_event(self.probe_move)
				self.img_widget.value = convert_to_bytes(self.img_norm, self.mag_slider.value, mode=self.tabs.selected_index, xy=self.coords, color=self.line_color_button.value)
			elif self.tabs.selected_index == 2: # Meas tab
				original_offset = [int(x[:-2]) for x in self.img_widget.layout.object_position.split()]
				self.coords = [(int((event['relativeX'] - original_offset[0])/self.mag_slider.value),
						int((event['relativeY'] - original_offset[1])/self.mag_slider.value))]
				self.mouse_move.on_dom_event(self.meas_move)
		if event['event']=='mouseup':
			# Disconnect proper events on releasing a mouse button
			if self.tabs.selected_index == 0: # App tab
				self.start_coord = None
				self.mouse_move.on_dom_event(self.pan, remove=True)
			elif self.tabs.selected_index == 1: # Probe tab
				self.mouse_move.on_dom_event(self.probe_move, remove=True)
			elif self.tabs.selected_index == 2: # Meas tab
				self.mouse_move.on_dom_event(self.meas_move, remove=True)

	def set_probe_size(self, change):
		"""Changes the size of the rectangle self.probe based on the slider self.probe_size_widget"""
		if change['new'] < 1:
			self.probe_size_widget.value = 1
		# Probe size only effects displays when we click again to redraw probe so nothing else is needed here

	@output.capture()
	def probe_move(self, event):
		"""Moves the probe, redraws the image, and updates the probe tab on a move event. It is throttled"""
		original_offset = [int(x[:-2]) for x in self.img_widget.layout.object_position.split()]
		center = (int((event['relativeX'] - original_offset[0])/self.mag_slider.value),
				int((event['relativeY'] - original_offset[1])/self.mag_slider.value))
		self.coords = [(center[0] - (self.probe_size_widget.value // 2), center[1] - (self.probe_size_widget.value // 2)),
				(center[0] + (self.probe_size_widget.value // 2), center[1] + (self.probe_size_widget.value // 2))]
		self.img_widget.value = convert_to_bytes(self.img_norm, self.mag_slider.value, mode=self.tabs.selected_index, xy=self.coords, color=self.line_color_button.value)
		self.update_probe_tab(center)

	def update_probe_tab(self, center):
		"""Display the proper statistics inside the probe tab for where self.probe is in the image"""
		size_x,size_y = self.eman_object[self.index.value]['nx'], self.eman_object[self.index.value]['ny']
		# According to coords/center 0,0 is top left but we want it to be bottom left
		x_corner = self.coords[0][0]
		y_corner = size_y - self.coords[0][1]
		x_center = center[0]
		y_center = size_y - center[1]
		clipped_img = self.eman_object[self.index.value].get_clip(Region(x_corner, y_corner, self.probe_size_widget.value, self.probe_size_widget.value))
		if x_center >= size_x or x_center < 0:
			val = 0
		elif y_center >= size_y or y_center < 0:
			val = 0
		else:
			val = self.cur_image[x_center, y_center]
		left_probe_col = widget.VBox([widget.Label(value=f"Point Value (ctr pix): {val:1.3f}"),
								widget.Label(value=f"Area Avg: {clipped_img['mean']:1.3f}"),
								widget.Label(value=f"Area Sig: {clipped_img['sigma']:1.3f}"),
								widget.Label(value=f"Skewness: {clipped_img['skewness']:1.3f}"),
								widget.Label(value=f"Kurtosis: {clipped_img['kurtosis']:1.3f}")],
								layout=text_layout())
		right_probe_col = widget.VBox([self.line_color_button,
								widget.Label(value=f"Area Avg (!=0): {clipped_img['mean_nonzero']:1.3f}"),
								widget.Label(value=f"Area Sig (!=0): {clipped_img['sigma_nonzero']:1.3f}"),
								widget.Label(value=f"Center Coord: {x_center}, {y_center}"),
								widget.Label(value=f"dcen ({x_center-size_x//2}, {y_center-size_y//2})")],
								layout=text_layout())
		probe_tab = widget.VBox([self.probe_size_widget,
					widget.HBox([left_probe_col, right_probe_col])],
					layout=text_layout())
		self.tabs.children = [self.tabs.children[0], probe_tab, self.tabs.children[2]]

	def meas_move(self, event):
		"""Moves the endpoint of the measure line, redraws the image, and updates the Meas tab on a move event. It is throttled"""
		original_offset = [int(x[:-2]) for x in self.img_widget.layout.object_position.split()]
		self.coords = [self.coords[0],
				(int((event['relativeX'] - original_offset[0])/self.mag_slider.value),
				int((event['relativeY'] - original_offset[1])/self.mag_slider.value))]
		self.img_widget.value = convert_to_bytes(self.img_norm, self.mag_slider.value, mode=self.tabs.selected_index, xy=self.coords, color=self.line_color_button.value)
		self.update_meas_tab()

	def update_meas_tab(self):
		"""Display the proper statistics inside the meas tab for the line at self.coords in the image and the apix given by
		self.aPix_slider"""
		size_x,size_y = self.eman_object[self.index.value]['nx'], self.eman_object[self.index.value]['ny']
		# Subtracting flips the y as required by coords 0,0 is top left but we want it to be bottom left
		meas_start = (self.coords[0][0], size_y - self.coords[0][1])
		meas_end = (self.coords[1][0], size_y - self.coords[1][1])
		if self.fourier_buttons.value == "Real":
			dx = (meas_end[0] - meas_start[0])*self.aPix_slider.value
			dy = (meas_end[1] - meas_start[1])*self.aPix_slider.value
			if meas_end[0] >= size_x or meas_end[0] < 0:
				val = 0
			elif meas_end[1] >= size_y or meas_end[1] < 0:
				val = 0
			else:
				val = self.cur_image[meas_end[0], meas_end[1]]
			left_meas_col = widget.VBox([widget.Label(value=f"Start: {meas_start[0]},{meas_start[1]}"),
							widget.Label(value=f"dx,dx: {dx:1.2f} A, {dy:1.2f} A"),
							widget.Label(value=f"Value: {val:1.4g}")],
							layout=text_layout())
			right_meas_col = widget.VBox([widget.Label(value=f"End: {meas_end[0]},{meas_end[1]}"),
							widget.Label(value=f"Len: {hypot(dx, dy):1.3f} A")],
							layout=text_layout())
			meas_tab = widget.VBox([widget.HBox([self.aPix_slider, self.line_color_button]),
							widget.HBox([left_meas_col, right_meas_col], layout=text_layout())])
		else:
			dx = meas_end[0] - meas_start[0]
			dy = meas_end[1] - meas_start[1]
			p = [meas_start[0]-size_x//2, meas_start[1]-size_y//2, meas_end[0]-size_x//2, meas_end[1]-size_y//2]
			p = [round(x) for x in p]
			q = [p[0]/size_x, p[1]/size_y, p[2]/size_x, p[3]/size_y]
			q = [f"{(self.aPix_slider.value/x):.1f}A" if x !=0 else "INF" for x in q]
			w = [hypot(p[0]/size_x, p[1]/size_y), hypot(p[2]/size_x, p[3]/size_y)]
			w = [f"{(self.aPix_slider.value/x):.1f}A" if x !=0 else "INF" for x in w]
			size_fx, size_fy = self.fft.get_xsize(), self.fft.get_ysize()
			x, y = meas_end[0]+1, meas_end[1]
			if x < (size_fx // 2):
				x = (size_fx // 2) - x
			else:
				x -= (size_fx // 2)
			if y < (size_fy // 2):
				y += (size_fy // 2)
			else:
				y -= (size_fy // 2)

			if x >= size_x or x < 0:
				val = 0
			elif y >= size_y or y < 0:
				val = 0
			else:
				val = self.fft[x,y]
			left_meas_col = widget.VBox([widget.Label(value=f"Start: {p[0]},{p[1]} ({q[0]}, {q[1]})"),
							widget.Label(value=f"dx,dx: {dx:1.0f} px, {dy:1.0f} px")],
							layout=text_layout())
			right_meas_col = widget.VBox([widget.Label(value=f"End: {p[2]}, {p[3]} ({q[2]}, {q[3]})"),
							widget.Label(value=f"Len: {hypot(dx, dy):1.1f} px ({w[0]}-> {w[1]})")],
							layout=text_layout())
			meas_tab = widget.VBox([self.aPix_slider,
							widget.HBox([left_meas_col, right_meas_col], layout=text_layout()),
							widget.Label(value=f"Value: {val.real:1.4g} + {val.imag:1.4g} i @({x},{y})"),
							widget.Label(value=f"	({abs(val):1.4g} + {atan2(val.imag, val.real)*57.295779513:1.4g})")],
							layout=text_layout())
		self.tabs.children = [self.tabs.children[0], self.tabs.children[1], meas_tab]

	def display_control_panel(self, event):
		"""Toggles the control panel visibility on left click or middle click inside the plot"""
		self.controls_above = not self.controls_above #flip the bool
		self.output.clear_output()
		with self.output:
			if self.controls_above:
				contrast_box = widget.VBox([self.fourier_buttons, self.invert_button, self.auto_button, self.full_button],
								layout=widget.Layout(width='10%'))
				self.tabs.layout = widget.Layout(width='35%')
				self.fourier_buttons.style={'button_width': '99%'}
				# sliders = widget.VBox([self.mag_slider, self.min_slider, self.max_slider, self.brightness,
				#						self.contrast, self.index], layout=widget.Layout(width='33%'))
				sliders = widget.VBox([self.mag_slider, self.slider_tabs, self.index], layout=widget.Layout(width='33%'))
				controls = widget.HBox([self.tabs, contrast_box, sliders, self.hist_fig.canvas])
				controls.layout = controls_layout(self.controls_above)
				self.img_widget.layout = image_layout(self.controls_above)
				display(widget.VBox([controls, self.img_widget]))
			else:
				contrast_box = widget.HBox([self.invert_button, self.auto_button, self.full_button])
				# controls = widget.VBox([self.tabs, self.fourier_buttons, contrast_box, self.mag_slider,
				#						 self.min_slider, self.max_slider, self.brightness, self.contrast,
				#						 self.hist_fig.canvas, self.index])
				self.tabs.layout = widget.Layout(width='100%')
				self.slider_tabs.layout = widget.Layout(width='100%')
				self.fourier_buttons.style={'button_width': '32%'}
				controls = widget.VBox([self.tabs, self.fourier_buttons, contrast_box, self.mag_slider,
								self.slider_tabs, self.hist_fig.canvas, self.index])
				controls.layout = controls_layout(self.controls_above)
				self.img_widget.layout = image_layout(self.controls_above)
				display(widget.HBox([controls, self.img_widget]))
		display(self.output)
