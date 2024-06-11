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
import matplotlib.widgets as mplwidget
import matplotlib.patches as patches
from past.utils import old_div
import ipywidgets as widget
from IPython.display import display
import asyncio
from time import time
from matplotlib.backend_bases import MouseEvent, MouseButton


def controls_layout():
    """Sets the layout for the controls box"""
    return widget.Layout(width='36%', border='solid 1px black', margin='0px 10px 10px 0px', padding='5px 5px 5px 5px')

def image_layout():
    """Sets the layout for the image box"""
    # return widget.Layout(width='64%', border='solid 1px black', margin='0px 10px 10px 0px', padding='5px 5px 5px 5px')
    return widget.Layout(margin='0px 10px 10px 0px', padding='5px 5px 5px 5px')

def button_layout():
    """Sets the layout for buttons"""
    return widget.Layout(width='32%')

def slider_layout():
    """Sets the layout for sliders"""
    return widget.Layout(width="99%")

def text_layout():
    """Sets the layout for text"""
    return widget.Layout(justify_content='space-around', align_content='flex-start')


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

    def __init__(self, eman_img_stack, show_controls=False):
        """Initialize the object and create the plot and widgets
        Inputs:
            self
            eman_img_stack - a list of EmanData objects, one that "Show2d" would be an option in e2display.py
        Outputs:
            none, displays widget when called"""

        style={'description_width': 'initial', 'handle_color': 'grey'}

        # Define widgets
        self.eman_object = eman_img_stack
        self.index = widget.IntSlider(min=0, max=len(self.eman_object),value=0, description='N#', layout=slider_layout(), style=style)
        self.cur_image = self.eman_object[self.index.value].numpy()
        self.fft = None
        self.fft_display = self.fft
        self.image_min = self.eman_object[self.index.value]['minimum']
        self.image_max = self.eman_object[self.index.value]['maximum']
        self.mean = self.eman_object[self.index.value]['mean']
        self.sigma = self.eman_object[self.index.value]['sigma']
        self.min_slider = widget.FloatSlider(min=self.image_min, max=self.image_max, step=0.001, value=self.image_min, description="Min", readout_format=".3f", layout=slider_layout(), style=style)
        self.max_slider = widget.FloatSlider(min=self.image_min, max=self.image_max, step=0.001, value=self.image_max, description="Max", readout_format=".3f", layout=slider_layout(), style=style)
        self.brightness = widget.FloatSlider(min=-1, max=1, step=0.001, value=0, description="Brt", readout_format=".3f", layout=slider_layout(), style=style)
        self.contrast = widget.FloatSlider(min=0, max=1, step=0.001, value=0.5, description="Cont", readout_format=".3f", layout=slider_layout(), style=style)
        self.auto_button = widget.Button(description="autoC", layout=button_layout())
        self.full_button = widget.Button(description="fullC", layout=button_layout())
        self.invert_button = widget.ToggleButton(description="Invert", layout=button_layout())
        self.changing_min_max = False
        self.changing_brt_cont = False
        self.fourier_buttons = widget.ToggleButtons(options=['Real', 'Amp', 'Pha'], style={'button_width': '32%'})
        self.mag_slider = widget.FloatSlider(min=0.1, max=10.0, value=1, description="Mag", layout=slider_layout(), style=style)
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
        app_tab = widget.HBox([size, stats, apix],  layout=text_layout())
        # Probe tab
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
                                 widget.Label(value= "Point Value (ctr pix): "), 
                                 widget.HBox([left_probe_col, right_probe_col], layout=text_layout())])
        # Measure tab
        self.aPix_slider = widget.FloatSlider(value=self.apix_x , min=0.5, max=10, step=0.001, description="A/Pix", readout_format=".3f", layout=slider_layout(), style=style)
        left_meas_col = widget.VBox([widget.Label(value="Start: 0,0"),
                                     widget.Label(value="dx,dx: 0"),
                                     widget.Label(value="Value: ?")],
                                   layout=text_layout())
        right_meas_col = widget.VBox([widget.Label(value="End 0,0"),
                                      widget.Label(value="Len: 0")],
                                    layout=text_layout())
        meas_tab = widget.VBox([self.aPix_slider, widget.HBox([left_meas_col, right_meas_col], layout=text_layout())])
        self.tabs = widget.Tab()
        self.tabs.children = [app_tab, probe_tab, meas_tab]
        self.tabs.titles = ["App", "Probe", "Meas"]
        self.output = widget.Output()

        self.background = None
        self.meas_start = None
        self.creating_bg = False
        self.control_panel = show_controls


        # Create figure we will modify
        with plt.ioff():
            self.fig, self.axs = plt.subplots(1, 1)
            self.imgplt = self.axs.imshow(self.cur_image, cmap="Greys_r", origin='lower', interpolation='none', resample=False)
        plt.axis("off")
        plt.tight_layout(pad=0.01)
        self.fig.canvas.toolbar_position = "bottom"
        self.probe = self.axs.add_patch(patches.Rectangle((0,0), self.probe_size_widget.value, self.probe_size_widget.value, animated=True, fill=False, visible=False, 
                                                    edgecolor="olive", linestyle='-', linewidth=2))
        self.line , = self.axs.plot([0,1], [0,1], animated=True, visible=False, linestyle='-', linewidth=1.5, color='m')
        self.set_mode({'new':self.tabs.selected_index})
        # grab the background on every draw
        self.draw_event = self.fig.canvas.mpl_connect("draw_event", self.on_draw)
        self.left_click = self.fig.canvas.mpl_connect("button_press_event", self.display_control_panel)

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
        self.tabs.observe(self.set_mode, 'selected_index')
        self.probe_size_widget.observe(self.set_probe_size, 'value')
        self.mag_slider.observe(self.update_zoom, 'value')

        # Start with auto contrast
        self.auto_button.click()

        # Establish Layout
        with self.output:
            if self.control_panel:
                contrast_box = widget.HBox([self.invert_button, self.auto_button, self.full_button])
                controls = widget.VBox([self.tabs, self.fourier_buttons, contrast_box, self.mag_slider, 
                                        self.min_slider, self.max_slider, self.brightness, self.contrast,
                                        self.hist_fig.canvas, self.index])
                controls.layout = controls_layout()
                self.fig.canvas.layout = image_layout()
                display(widget.HBox([controls, self.fig.canvas]))
            else:
                display(self.fig.canvas)
        display(self.output)




    def __del__(self):
        """When object is deleted tell jupyter to close the figures"""
        plt.close(self.fig)
        plt.close(self.hist)


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
        app_tab = widget.HBox([size, stats, apix],  layout=text_layout())
        if self.tabs.selected_index > 0:
            self.tabs.selected_index = None
        self.tabs.children = [app_tab, self.tabs.children[1], self.tabs.children[2]]
        self.tabs.selected_index = selected_tab

    def change_image(self, change):
        """Updates the image displayed from the stack of images self.eman_object based on 
        the slider self.index. Will keep the previous fourier status. Also updates histogram"""
        self.cur_image = self.eman_object[change['new']].numpy()
        self.fft = None
        self.fft_display = None
        if self.fourier_buttons.value == 'Real':
            self.imgplt.set_data(self.cur_image)
            self.image_min = self.eman_object[change['new']]['minimum']
            self.image_max = self.eman_object[change['new']]['maximum']
            self.mean = self.eman_object[change['new']]['mean']
            self.sigma = self.eman_object[change['new']]['sigma']
            self.hist_axs.cla()
            self.hist_axs.hist(self.cur_image.flatten(), bins='auto')
        elif self.fourier_buttons.value == 'Amp':
            self.fft = self.eman_object[change['new']].do_fft()
            self.fft_display = self.fft.process("xform.fourierorigin.tocenter")
            self.fft_display = self.fft_display.get_fft_amplitude()
            self.image_min = self.fft_display['minimum']
            self.image_max = self.fft_display['maximum']
            self.mean = self.fft_display['mean']
            self.sigma = self.fft_display['sigma']
            self.fft_display = self.fft_display.numpy()
            self.imgplt.set_data(self.fft_display)
            self.hist_axs.cla()
            self.hist_axs.hist(self.fft_display.flatten(), bins='auto', range=(self.image_min, self.image_max))
        elif self.fourier_buttons.value == 'Pha':
            self.fft = self.eman_object[change['new']].do_fft()
            self.fft_display = self.fft.process("xform.fourierorigin.tocenter")
            self.fft_display = self.fft_display.get_fft_phase()
            self.image_min = self.fft_display['minimum']
            self.image_max = self.fft_display['maximum']
            self.mean = self.fft_display['mean']
            self.sigma = self.fft_display['sigma']
            self.fft_display = self.fft_display.numpy()
            self.imgplt.set_data(self.fft_display)
            self.hist_axs.cla()
            self.hist_axs.hist(self.fft_display.flatten(), bins='auto', range=(self.image_min, self.image_max))
        self.fig.canvas.draw_idle()

    def fourier_transform(self, change):
        """Switches between displaying the real image, the fourier amplitude, or the fourier phase of the image
        based on the toggle buttons self.fourier_buttons. Also changes histogram."""
        if change['new'] == 'Real':
            self.imgplt.set_data(self.cur_image)
            self.image_min = self.eman_object[self.index.value]['minimum']
            self.image_max = self.eman_object[self.index.value]['maximum']
            self.mean = self.eman_object[self.index.value]['mean']
            self.sigma = self.eman_object[self.index.value]['sigma']
            self.hist_axs.cla()
            self.hist_axs.hist(self.cur_image.flatten(), bins='auto')
        elif change['new'] == 'Amp':
            if self.fft == None:
                self.fft = self.eman_object[self.index.value].do_fft()
            self.fft_display = self.fft.process("xform.fourierorigin.tocenter")
            self.fft_display = self.fft_display.get_fft_amplitude()
            self.image_min = self.fft_display['minimum']
            self.image_max = self.fft_display['maximum']
            self.mean = self.fft_display['mean']
            self.sigma = self.fft_display['sigma']
            self.fft_display = self.fft_display.numpy()
            self.imgplt.set_data(self.fft_display)
            self.hist_axs.cla()
            self.hist_axs.hist(self.fft_display.flatten(), bins='auto', range=(self.image_min, self.image_max))
        elif change['new'] == 'Pha':
            if self.fft == None:
                self.fft = self.eman_object[self.index.value].do_fft()
            self.fft_display = self.fft.process("xform.fourierorigin.tocenter")
            self.fft_display = self.fft_display.get_fft_phase()
            self.image_min = self.fft_display['minimum']
            self.image_max = self.fft_display['maximum']
            self.mean = self.fft_display['mean']
            self.sigma = self.fft_display['sigma']
            self.fft_display = self.fft_display.numpy()
            self.imgplt.set_data(self.fft_display)
            self.hist_axs.cla()
            self.hist_axs.hist(self.fft_display.flatten(), bins='auto', range=(self.image_min, self.image_max))
        self.auto_button.click()
        self.fig.canvas.draw_idle()

    def autoC(self, b):
        """Sets the self.min_slider and self.max_slider values for auto contrast"""
        if self.fourier_buttons.value == 'Real':
            self.min_slider.min = self.image_min
            self.min_slider.max = self.image_max
            self.max_slider.min = self.image_min
            self.max_slider.max = self.image_max
            self.min_slider.value = max(self.image_min, self.mean-3.0*self.sigma)
            self.max_slider.value = min(self.image_max, self.mean+3.0*self.sigma)
        else:
            self.min_slider.min = 0
            self.min_slider.max = min(self.image_max, self.mean+20.0*self.sigma)
            self.max_slider.min = 0
            self.max_slider.max = min(self.image_max, self.mean+20.0*self.sigma)
            self.min_slider.value = 0
            self.max_slider.value = min(self.image_max, self.mean+4.0*self.sigma)

    def fullC(self, b):
        """Sets the self.min_slider and self.max_slider values for auto contrast"""
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
            self.min_slider.max = min(self.image_max, self.mean+20.0*self.sigma)
            self.max_slider.max = min(self.image_max, self.mean+20.0*self.sigma)
            self.min_slider.value = 0
            self.max_slider.value = self.image_max

    def invert(self, change):
        """Inverts the image"""
        if change['new']:
            self.imgplt.set_cmap("Greys")
        else:
            self.imgplt.set_cmap("Greys_r")
        self.fig.canvas.draw_idle()

    @throttle(0.3)
    def update_brt_cont(self, min, max):
        """Updates self.brightness and self.contrast sliders to reflect changes in self.min_slider or
        self.max_slider. Throttled to prevent jumpiness."""
        self.brightness.value = -0.5*(min+max-(self.image_max+self.image_min))/((self.image_max-self.image_min))
        self.contrast.value = 1.0 - old_div((min-max), (2.0*(self.image_min-self.image_max)))

    @throttle(0.3)
    def update_min_max(self, brt, cont):
        """Updates self.min_slider and self.max_slider to reflect changes in self.brightness or
        self.contrast sliders. Throttled to prevent jumpiness."""
        m0 = old_div((self.image_min+self.image_max),2.0)-(self.image_max-self.image_min)*(1.0-cont)-brt*(self.image_max-self.image_min)
        m1 = old_div((self.image_min+self.image_max),2.0)+(self.image_max-self.image_min)*(1.0-cont)-brt*(self.image_max-self.image_min)
        # with self.min_slider.hold_sync():
        self.min_slider.min = min(m0, self.min_slider.min)
        # with self.max_slider.hold_sync():
        self.max_slider.max = max(m1, self.max_slider.max)
        self.min_slider.value = m0
        self.max_slider.value = m1

    def update_min(self, change):
        """Updates the minimum value mapped to the colormap based on a change to self.min_slider.
        If the change was due to changing min_slider itself and not from changing self.brightness or self.contrast,
        also updates those sliders."""
        if not self.changing_min_max:
            # Update Brightness and contrast
            self.changing_min_max = True
            self.update_brt_cont(change['new'], self.max_slider.value)
            self.changing_min_max = False
        # Update the image's colormap
        self.imgplt.norm.vmin = change['new']
        # Redraw the figure to ensure it updates
        self.fig.canvas.draw_idle()

    def update_max(self, change):
        """Updates the maximum value mapped to the colormap based on a change to self.max_slider.
        If the change was due to changing max_slider itself and not from changing self.brightness or self.contrast,
        also updates those sliders."""
        if not self.changing_min_max:
            # Update Brightness and contrast
            self.changing_min_max = True
            self.update_brt_cont(self.min_slider.value, change['new'])
            self.changing_min_max = False
        # Update the image's colormap
        self.imgplt.norm.vmax = change['new']
        # Redraw the figure to ensure it updates
        self.fig.canvas.draw_idle()

    def update_brt(self, change):
        """Updates self.min_slider and self.max_slider when the user changes self.brightness."""
        if not self.changing_brt_cont:
            self.changing_brt_cont = True
            self.update_min_max(change['new'], self.contrast.value)
            self.changing_brt_cont = False
        self.fig.canvas.draw_idle()

    def update_cont(self, change):
        """Updates self.min_slider and self.max_slider when the user changes self.contrast."""
        if not self.changing_brt_cont:
            self.changing_brt_cont = True
            self.update_min_max(self.brightness.value, change['new'])
            self.changing_brt_cont = False
        self.fig.canvas.draw_idle()

    def update_zoom(self, change):
        """Updates the image based on a change in self.mag_slider. A 2x increase in the value will be a 2x increase
        in the image side length"""
        scale_factor = change['old']/change['new']
        cur_xlim = self.axs.get_xlim()
        cur_ylim = self.axs.get_ylim()
        # xdata = old_div(cur_xlim[1]-cur_xlim[0], 2)
        # ydata = old_div(cur_ylim[1]-cur_ylim[0], 2)
        xdata = old_div(self.eman_object[self.index.value].get_xsize(), 2)
        ydata = old_div(self.eman_object[self.index.value].get_ysize(), 2)
        x_left = xdata - cur_xlim[0]
        x_right = cur_xlim[1] - xdata
        y_top = ydata - cur_ylim[0]
        y_bottom = cur_ylim[1] - ydata
        self.axs.set_xlim([xdata - x_left*scale_factor,
                           xdata + x_right*scale_factor])
        self.axs.set_ylim([ydata - y_top*scale_factor,
                           ydata + y_bottom*scale_factor])
        self.fig.canvas.draw_idle()

    def set_probe_size(self, change):
        """Changes the size of the rectangle self.probe based on the slider self.probe_size_widget"""
        if change['new'] < 1:
            self.probe_size_widget.value = 1
        self.probe.set(height=self.probe_size_widget.value, width=self.probe_size_widget.value)

    def set_mode(self, change):
        """Based off of which tab is selected (App, Probe, or Meas), sets self.probe and self.line
        to the correct visibility and connects button_press_event and button_release_event to the 
        correct functions."""
        if change['new'] == 0:
            self.probe.set_visible(False)
            self.line.set_visible(False)
            try:
                self.fig.canvas.mpl_disconnect(self.probe_click)
                self.fig.canvas.mpl_disconnect(self.probe_release)
            except: pass
            try:
                self.fig.canvas.mpl_disconnect(self.meas_click)
                self.fig.canvas.mpl_disconnect(self.meas_release)
            except: pass
        if change['new'] == 1:
            self.probe.set_visible(True)
            self.line.set_visible(False)
            self.probe_click = self.fig.canvas.mpl_connect('button_press_event', self.probe_mouse_press)
            #self.probe_move = self.canvas.mpl_connect('motion_notify_event', self.probe_mouse_move)
            self.probe_release = self.fig.canvas.mpl_connect('button_release_event', self.probe_mouse_release)
            try:
                self.fig.canvas.mpl_disconnect(self.meas_click)
                self.fig.canvas.mpl_disconnect(self.meas_release)
            except: pass
        if change['new'] == 2:
            self.probe.set_visible(False)
            self.line.set_visible(True)
            self.meas_click = self.fig.canvas.mpl_connect('button_press_event', self.meas_mouse_press)
            self.meas_release = self.fig.canvas.mpl_connect('button_release_event', self.meas_mouse_release)
            try:
                self.fig.canvas.mpl_disconnect(self.probe_click)
                self.fig.canvas.mpl_disconnect(self.probe_release)
            except: pass

    # Functions for matplotlib plot interaction
    def display_control_panel(self, event):
        """Toggles the control panel visibility on left click inside the plot"""
        if self.fig.canvas.manager.toolbar.mode != '':
            return
        if event.inaxes != self.axs:
            return
        if event.button is MouseButton.RIGHT:
            self.control_panel = not self.control_panel #flip the bool
            self.output.clear_output()
            self.tabs.selected_index = 0
            with self.output:
                if self.control_panel:
                    contrast_box = widget.HBox([self.invert_button, self.auto_button, self.full_button])
                    controls = widget.VBox([self.tabs, self.fourier_buttons, contrast_box, self.mag_slider, 
                                            self.min_slider, self.max_slider, self.brightness, self.contrast,
                                            self.hist_fig.canvas, self.index])
                    controls.layout = controls_layout()
                    self.fig.canvas.layout = image_layout()
                    display(widget.HBox([controls, self.fig.canvas]))
                else:
                    display(self.fig.canvas)
            # change display to show control panel

    def on_draw(self, event):
        """When a draw event occurs, re-capture the background without self.probe or self.line"""
        if self.creating_bg:
            return
        self.creating_bg = True
        self.probe.set_visible(False)
        self.line.set_visible(False)
        self.fig.canvas.draw_idle()
        self.background = self.fig.canvas.copy_from_bbox(self.fig.bbox)
        if self.tabs.selected_index == 1:
            self.probe.set_visible(True)
            self.draw_probe()
        if self.tabs.selected_index == 2:
            self.line.set_visible(True)
            self.draw_meas()
        self.creating_bg = False

    def draw_probe(self):
        """Draw self.probe using blitting on top of the canvas"""
        self.fig.canvas.restore_region(self.background)
        self.fig.draw_artist(self.probe)
        self.fig.canvas.blit(self.fig.bbox)
        # self.canvas.flush_events()

    @throttle(0.4)
    def draw_probe_throttled(self):
        """A throttled version of draw_probe used when drawing off of mouse_move events because blitting is not
        supported with the ipympl backend and we get afterimages. Reducing the calls reduces the chance of getting an
        afterimage but I have still seen it happen"""
        self.fig.canvas.restore_region(self.background)
        self.fig.draw_artist(self.probe)
        self.fig.canvas.blit(self.fig.bbox)
        # self.canvas.flush_events()

    def probe_mouse_press(self, event):
        """On a left mouse press inside the axis, move the probe to that location and connect mouse_move events
        for the probe. This only occurs when no other mode is selected in the matplotlib toolbar."""
        # Trigger probe rectangle to move to cursor and follow so long as button is clicked
        if self.fig.canvas.manager.toolbar.mode != '':
            return
        if event.inaxes != self.axs:
            return
        if event.button is MouseButton.LEFT:
            # data coords: event.xdata, event.ydata
            # pixel coords event.x, event.y
            self.probe_move = self.fig.canvas.mpl_connect('motion_notify_event', self.probe_mouse_move)
            self.probe.set_x(event.xdata-old_div(self.probe_size_widget.value, 2))
            self.probe.set_y(event.ydata-old_div(self.probe_size_widget.value, 2))
            self.draw_probe()

    def probe_mouse_move(self, event):
        """When the mouse moves while we are using probe mode, move the rectangle to the mouse, re-draw the probe,
        and update the probe tab with statistics"""
        if event.inaxes != self.axs:
            return
        self.probe.set_x(event.xdata-old_div(self.probe_size_widget.value, 2))
        self.probe.set_y(event.ydata-old_div(self.probe_size_widget.value, 2))
        self.draw_probe_throttled()
        self.update_probe_tab()

    def probe_mouse_release(self, event):
        """When the left mouse button is released disconnect the move event and draw the probe one last time.
        This only occurs when no other mode is selected in the matplotlib toolbar."""
        if self.fig.canvas.manager.toolbar.mode != '':
            return
        if event.button is MouseButton.LEFT:
            self.draw_probe()
            self.fig.canvas.mpl_disconnect(self.probe_move)

    def update_probe_tab(self):
        """Display the proper statistics inside the probe tab for where self.probe is in the image"""
        # selected_tab = self.tabs.selected_index
        x_corner = int(round(self.probe.get_x()))
        y_corner = int(round(self.probe.get_y()))
        x_center, y_center = self.probe.get_center()
        x_center = int(round(x_center))
        y_center = int(round(y_center))
        size_x,size_y = self.eman_object[self.index.value]['nx'], self.eman_object[self.index.value]['ny']
        clipped_img = self.eman_object[self.index.value].get_clip(Region(x_corner, y_corner, self.probe_size_widget.value, self.probe_size_widget.value))
        left_probe_col = widget.VBox([widget.Label(value=f"Area Avg: {clipped_img['mean']:1.3f}"),
                                      widget.Label(value=f"Area Sig: {clipped_img['sigma']:1.3f}"),
                                      widget.Label(value=f"Skewness: {clipped_img['skewness']:1.3f}"),
                                      widget.Label(value=f"Kurtosis: {clipped_img['kurtosis']:1.3f}")],
                                    layout=text_layout())
        right_probe_col = widget.VBox([widget.Label(value=f"Area Avg (!=0): {clipped_img['mean_nonzero']:1.3f}"),
                                       widget.Label(value=f"Area Sig (!=0): {clipped_img['sigma_nonzero']:1.3f}"),
                                       widget.Label(value=f"Center Coord: {x_center}, {y_center}"),
                                       widget.Label(value=f"dcen ({x_center-size_x//2}, {y_center-size_y//2})")],
                                     layout=text_layout())
        val = self.cur_image[x_center, y_center]
        probe_tab = widget.VBox([self.probe_size_widget, 
                                 widget.Label(value=f"Point Value (ctr pix): {val:1.3f}"), 
                                 widget.HBox([left_probe_col, right_probe_col], layout=text_layout())])
        self.tabs.children = [self.tabs.children[0], probe_tab, self.tabs.children[2]]
        # self.tabs.selected_index = selected_tab

    def draw_meas(self):
        """Draw self.line using blitting on top of the canvas"""
        self.fig.canvas.restore_region(self.background)
        self.axs.draw_artist(self.line)
        self.fig.canvas.blit()
        # self.canvas.flush_events()

    @throttle(0.4)
    def draw_meas_throttled(self):
        """A throttled version of draw_meas used when drawing off of mouse_move events because blitting is not
        supported with the ipympl backend and we get afterimages. Reducing the calls reduces the chance of getting an
        afterimage but I have still seen it happen"""
        self.fig.canvas.restore_region(self.background)
        self.axs.draw_artist(self.line)
        self.fig.canvas.blit()

    def meas_mouse_press(self, event):
        """On a left mouse press inside the axis, save the initial point as the starting point for self.line
        and connect mouse_move events for measuring. This only occurs when no other mode is selected in the matplotlib toolbar."""
        if self.fig.canvas.manager.toolbar.mode != '':
            return
        if event.inaxes != self.axs:
            return
        if event.button is MouseButton.LEFT:
            self.meas_move = self.fig.canvas.mpl_connect('motion_notify_event', self.meas_mouse_move)
            self.meas_start = (event.xdata, event.ydata) # Is this right or should I be using pixel vals?

    def meas_mouse_move(self, event):
        """When the mouse moves while we are using meas mode, move the endpoint of the line to the mouse, re-draw the
        line using the saved starting point, and update the meas tab with statistics,"""
        if event.inaxes != self.axs:
            return
        self.line.set_data([self.meas_start[0], event.xdata], [self.meas_start[1], event.ydata])
        self.draw_meas_throttled()
        self.update_meas_tab((event.xdata, event.ydata))

    def meas_mouse_release(self, event):
        """When the left mouse button is released disconnect the move event and draw the line one last time.
        This only occurs when no other mode is selected in the matplotlib toolbar."""
        if self.fig.canvas.manager.toolbar.mode != '':
            return
        if event.inaxes != self.axs:
            return
        if event.button is MouseButton.LEFT:
            self.fig.canvas.mpl_disconnect(self.meas_move)
            self.line.set_data([self.meas_start[0], event.xdata], [self.meas_start[1], event.ydata])
            self.draw_meas()
            self.update_meas_tab((event.xdata, event.ydata))
            # self.event=event

    def update_meas_tab(self, meas_end):
        """Display the proper statistics inside the meas tab for where self.line is in the image and the apix given by 
        self.aPix_slider"""
        # selected_tab = self.tabs.selected_index
        if self.fourier_buttons.value == "Real":
            dx = (meas_end[0] - self.meas_start[0])*self.aPix_slider.value
            dy = (meas_end[1] - self.meas_start[1])*self.aPix_slider.value
            val = self.cur_image[int(round(meas_end[0])), int(round(meas_end[1]))]
            left_meas_col = widget.VBox([widget.Label(value=f"Start: {int(round(self.meas_start[0]))},{int(round(self.meas_start[1]))}"),
                                         widget.Label(value=f"dx,dx: {dx:1.2f} A, {dy:1.2f} A"),
                                         widget.Label(value=f"Value: {val:1.4g}")],
                                       layout=text_layout())
            right_meas_col = widget.VBox([widget.Label(value=f"End: {int(round(meas_end[0]))},{int(round(meas_end[1]))}"),
                                          widget.Label(value=f"Len: {hypot(dx, dy):1.3f} A")],
                                        layout=text_layout())
            meas_tab = widget.VBox([self.aPix_slider, widget.HBox([left_meas_col, right_meas_col], layout=text_layout())])
        else:
            dx = meas_end[0] - self.meas_start[0]
            dy = meas_end[1] - self.meas_start[1]
            size_x,size_y = self.eman_object[self.index.value]['nx'], self.eman_object[self.index.value]['ny']
            p = [self.meas_start[0]-size_x//2, self.meas_start[1]-size_y//2, meas_end[0]-size_x//2, meas_end[1]-size_y//2]
            p = [round(x) for x in p]
            q = [p[0]/size_x, p[1]/size_y, p[2]/size_x, p[3]/size_y]
            q = [f"{(self.aPix_slider.value/x):.1f}A" if x !=0 else "INF" for x in q]
            w = [hypot(p[0]/size_x, p[1]/size_y), hypot(p[2]/size_x, p[3]/size_y)]
            w = [f"{(self.aPix_slider.value/x):.1f}A" if x !=0 else "INF" for x in w]
            size_fx, size_fy = self.fft.get_xsize(), self.fft.get_ysize()
            x, y = int(meas_end[0])+1, int(meas_end[1])
            if x < old_div(size_fx, 2): x = old_div(size_fx, 2) - x
            else: x -= old_div(size_fx, 2)
            if y < old_div(size_fy, 2): y += old_div(size_fy, 2)
            else: y -= old_div(size_fy, 2)
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
                                    widget.Label(value=f"       ({abs(val):1.4g} + {atan2(val.imag, val.real)*57.295779513:1.4g})")],
                                  layout=text_layout())
        self.tabs.children = [self.tabs.children[0], self.tabs.children[1], meas_tab]
        # self.tabs.selected_index = selected_tab

