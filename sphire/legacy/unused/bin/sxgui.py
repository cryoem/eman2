self.customButtonStyle = """1
			SXLogoButton {{background-color: rgba(0, 0, 0, 0); border: 0px solid black; border-radius: 0px; image: url("{0}");}}
			SXLogoButton:focus {{background-color: rgba(0, 0, 0, 0); border: 0px solid grey; border-radius: 0px; image: url("{0}");}}
			SXLogoButton:pressed {{background-color: rgba(0, 0, 0, 0); border: 0px solid red; border-radius: 0px; image: url("{0}");}}
			""".format(logo_file_path)
self.customButtonStyleClicked = """2
			SXLogoButton {{background-color: rgba(0, 0, 0, 0); border: 0px solid black; border-radius: 0px; image: url("{0}");}}
			SXLogoButton:focus {{background-color: rgba(0, 0, 0, 0); border: 0px solid grey; border-radius: 0px; image: url("{0}");}}
			SXLogoButton:pressed {{background-color: rgba(0, 0, 0, 0); border: 0px solid red; border-radius: 0px; image: url("{0}");}}
			""".format(logo_file_path)
self.customButtonStyle = """3
			SXPictogramButton {{background-color: rgba(0, 0, 0, 0); border: 2px solid rgba(0, 0, 0, 0); border-radius: {1}px; image: url("{0}");}}
			SXPictogramButton:focus {{background-color: rgba(0, 0, 0, 0); border: 2px solid grey; border-radius: {1}px; image: url("{0}");}}
			SXPictogramButton:pressed {{background-color: rgba(0, 0, 0, 0); border: 2px solid rgb(153, 153, 153); border-radius: {1}px; image: url("{0}");}}
			""".format(pictogram_file_path, pictogram_width / 6)
self.customButtonStyleClicked = """4
			SXPictogramButton:pressed {{background-color: rgba(0, 0, 0, 0); border: 2px solid rgb(153, 153, 153); border-radius: {1}px; image: url("{0}");}}
			SXPictogramButton {{background-color: rgba(0, 0, 0, 0); border: 2px solid rgb(220, 220, 220); border-radius: {1}px; image: url("{0}");}}
			""".format(pictogram_file_path, pictogram_width / 6)
tab_widget.setStyleSheet("""QTabWidget::pane {5
			border-top: 2px solid #C2C7CB;
			position: absolute;
			top: -0.5em;
		}

		QTabWidget::tab-bar {
			alignment: center;
		}

		QTabBar::tab {
			background: qlineargradient(x1: 0, y1: 0, x2: 0, y2: 1,
					stop: 0 #E1E1E1, stop: 0.4 #DDDDDD,
					stop: 0.5 #D8D8D8, stop: 1.0 #D3D3D3);
			border: 2px solid #C4C4C3;
			border-bottom-color: #C2C7CB; /* same as the pane color */
			border-top-left-radius: 4px;
			border-top-right-radius: 4px;
			min-width: 8ex;
			padding: 2px;
		}

		QTabBar::tab:selected, QTabBar::tab:hover {
			background: qlineargradient(x1: 0, y1: 0, x2: 0, y2: 1,
					stop: 0 #fafafa, stop: 0.4 #f4f4f4,
					stop: 0.5 #e7e7e7, stop: 1.0 #fafafa);
		}

		QTabBar::tab:selected {
			border-color: #9B9B9B;
			border-bottom-color: #C2C7CB; /* same as pane color */
		}""")
"""6
#	def show_output_info(self):
#		QMessageBox.information(self, "sx* output","outdir is the name of the output folder specified by the user. If it does not exist, the directory will be created. If it does exist, the program will crash and an error message will come up. Please change the name of directory and restart the program.")
	"""
"""7
# ========================================================================================
class SXDialogCalculator(QDialog):
	def __init__(self, parent = None):
		super(QDialog, self).__init__(parent)
		
		self.setWindowModality(Qt.ApplicationModal)
		
		# self.setWindowTitle()
		self.setWindowTitle("Absolute Frequency Calculator")

		temp_label = QLabel("Calculate absolute frequency [1/Pixel] from resolution [A]", self)
		temp_label.move(50,50)
		
		# Create label widget
		temp_label = QLabel("Resolution [A]", self)
		# temp_label.setMinimumWidth(token_label_min_width)
		# grid_layout.addWidget(temp_label, grid_row, grid_col_origin, token_label_row_span, token_label_col_span)
		temp_label.move(50,100)
		self.edit_res = QLineEdit(self)
		self.edit_res.setText('Enter Resolution Here')
		self.edit_res.move(200,100)

		temp_label = QLabel("Pixel Size [A/Pixel]", self)
		# temp_label.setMinimumWidth(token_label_min_width)
		# grid_layout.addWidget(temp_label, grid_row, grid_col_origin, token_label_row_span, token_label_col_span)
		temp_label.move(50,200)
		self.edit_apix = QLineEdit(self)
		self.edit_apix.setText('Enter Pixel Size Here')
		self.edit_apix.move(200,200)
		
		self.btn_apply = QPushButton("Apply", self)
		self.btn_apply.move(50,300)
		self.btn_cancel = QPushButton("Cancel", self)
		self.btn_cancel.move(200,300)
		# self.connect(cmd_token_restore_widget[widget_index], SIGNAL("clicked()"), partial(self.handle_restore_widget_event, cmd_token, widget_index))
		
		### self.show()
"""
