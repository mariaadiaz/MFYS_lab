import numpy as np
from PyQt5 import QtWidgets

from vispy import app,scene
from vispy.scene import SceneCanvas, visuals
from vispy.app import use_app


CANVAS_SIZE = (300, 400)  # (width, height)
NUM_LINE_POINTS = 1000
H = np.zeros((4, 1000))
H2 = np.zeros((4, 1000))
ei = np.zeros([1, 100])
y_pred = np.zeros((100))
y_pred2 = np.zeros((100))
x_range = np.linspace(75,126,100)
pos = np.array([[75,0],[80,0],[100,0]])
y = np.zeros(len(x_range))
y1 = y + 0.5
y2 = y - 0.5
COLORMAP_CHOICES = ["viridis", "reds", "blues"]
LINE_COLOR_CHOICES = ["black", "blue", "green"]
n = 4

class CanvasWrapper():
    def __init__(self):
        app.use_app("PyQt5")
        self.canvas = scene.SceneCanvas(size=CANVAS_SIZE, keys='interactive')
        self.grid = self.canvas.central_widget.add_grid(spacing = 0)

        """Bayesian optimziation panel"""
        self.view_top1 = self.grid.add_view(0, 1, col_span = 2, bgcolor='white')
        self.view_top1.camera = "panzoom"
        self.view_top1.camera.set_range(x=(70, 135), y=(2, 4.8))
        sf_vals = pos

        line_pred = generate_optimization(x_range,y)
        line_pred2 = generate_optimization(x_range,y_pred)     
        vertices = fill_region(x_range,y1,y2)
        self.markers = visuals.Markers(pos = sf_vals, edge_color = 'black', face_color = 'black', size = 0.01, parent=self.view_top1.scene)
        self.themarker = visuals.Markers(pos = np.array([[0,0]]), edge_color = 'black', face_color = 'black', size = 0.01, parent=self.view_top1.scene)
        self.linepred = visuals.Line(line_pred, parent=self.view_top1.scene, color=LINE_COLOR_CHOICES[0], width = 0.01)
        self.linepred2 = visuals.Line(line_pred2, parent=self.view_top1.scene, color=LINE_COLOR_CHOICES[1], width = 0.01)
        self.area = visuals.Polygon(pos = vertices, color = 'lightblue', parent=self.view_top1.scene)

        yaxis = scene.AxisWidget(orientation='left', axis_label='Obj. func. (Model)', axis_font_size=8, axis_label_margin=30, tick_label_margin=2)
        yaxis.width_max = 50
        self.grid.add_widget(yaxis, row=0, col=0)
        yaxis.link_view(self.view_top1)

        x_axis = scene.AxisWidget(orientation='bottom', axis_label='Step frequencies', axis_font_size=8, axis_label_margin=20,tick_label_margin=15)
        x_axis.stretch = (1,0.2)
        self.grid.add_widget(x_axis, row=1, col=1, col_span=2)
        x_axis.link_view(self.view_top1)

        """Expected improvement panel"""
        self.view_top2 = self.grid.add_view(2, 1, col_span = 2, bgcolor='white')
        self.view_top2.camera = "panzoom"
        self.view_top2.camera.set_range(x=(70, 135), y=(-0.001, 0.1))

        line_ei = generate_optimization(x_range,ei)
        self.line_ei = visuals.Line(line_ei, parent=self.view_top2.scene, color='black', width=2)

        yaxis2 = scene.AxisWidget(orientation='left', axis_label='Acquisition func.', axis_font_size=8, axis_label_margin=30, tick_label_margin=2)
        yaxis2.width_max = 50
        self.grid.add_widget(yaxis2, row=2, col=0)
        yaxis2.link_view(self.view_top2)

        x_axis2 = scene.AxisWidget(orientation='bottom', axis_label='Step frequencies', axis_font_size=8, axis_label_margin=20,tick_label_margin=15)
        x_axis2.stretch = (1,0.2)
        self.grid.add_widget(x_axis2, row=3, col=1, col_span=2)
        x_axis2.link_view(self.view_top2)

        # Add gridlines for both previous panels#
        thegrid1= scene.visuals.GridLines(parent=self.view_top1.scene, color = 'black', scale= (0.01,0.01) )
        thegrid1.set_gl_state('translucent')
        thegrid2= scene.visuals.GridLines(parent=self.view_top2.scene, color = 'black', scale= (0.01,0.01) )
        thegrid2.set_gl_state('translucent')

        """PANEL 2: Coefficients of activation""" 

        self.view_bot = self.grid.add_view(4, 1, bgcolor='white')
        line_data = generate_compare_coefficients(NUM_LINE_POINTS,H[0])
        line_data2 = generate_compare_coefficients(NUM_LINE_POINTS,H2[0])
        self.line = visuals.Line(line_data, parent=self.view_bot.scene, color=LINE_COLOR_CHOICES[0], width = 2)
        self.line2 = visuals.Line(line_data2, parent=self.view_bot.scene, color=LINE_COLOR_CHOICES[1], width = 2)
        self.view_bot.camera = "panzoom"
        self.view_bot.camera.set_range(x=(0, 1), y=(0, 0.6))

        self.view_bot2 = self.grid.add_view(4, 2, bgcolor='white')
        line_data3 = generate_compare_coefficients(NUM_LINE_POINTS,H[1])
        line_data4 = generate_compare_coefficients(NUM_LINE_POINTS,H2[1])
        self.line3 = visuals.Line(line_data3, parent=self.view_bot2.scene, color=LINE_COLOR_CHOICES[0], width = 2)
        self.line4 = visuals.Line(line_data4, parent=self.view_bot2.scene, color=LINE_COLOR_CHOICES[1], width = 2)
        self.view_bot2.camera = "panzoom"
        self.view_bot2.camera.set_range(x=(0, 1), y=(0, 0.6))         # x and y are used to set the range of the axes

        yaxis3 = scene.AxisWidget(orientation='left', axis_label='Ci', axis_font_size=8, axis_label_margin=30, tick_label_margin=2)
        yaxis3.width_max = 50
        self.grid.add_widget(yaxis3, row=4, col=0)
        yaxis3.link_view(self.view_bot)
        x_axis3 = scene.AxisWidget(orientation='bottom', axis_label='Time (%)', axis_font_size=8, axis_label_margin=20,tick_label_margin=15)
        x_axis3.stretch = (1,0.2)
        self.grid.add_widget(x_axis3, row=5, col=1)
        x_axis3.link_view(self.view_bot)
        x_axis4 = scene.AxisWidget(orientation='bottom', axis_label='Time (%)', axis_font_size=8, axis_label_margin=20,tick_label_margin=15)
        x_axis4.stretch = (1,0.2)
        self.grid.add_widget(x_axis4, row=5, col=2)
        x_axis4.link_view(self.view_bot2)

        self.view_bot3 = self.grid.add_view(6, 1, bgcolor='white')
        line_data5 = generate_compare_coefficients(NUM_LINE_POINTS,H[2])
        line_data6 = generate_compare_coefficients(NUM_LINE_POINTS,H2[2])
        self.line5 = visuals.Line(line_data5, parent=self.view_bot3.scene, color=LINE_COLOR_CHOICES[0], width = 2)
        self.line6 = visuals.Line(line_data6, parent=self.view_bot3.scene, color=LINE_COLOR_CHOICES[1], width = 2)
        self.view_bot3.camera = "panzoom"
        self.view_bot3.camera.set_range(x=(0, 1), y=(0, 0.6))         # x and y are used to set the range of the axes

        self.view_bot4 = self.grid.add_view(6, 2, bgcolor='white')
        line_data7 = generate_compare_coefficients(NUM_LINE_POINTS,H[0])
        line_data8 = generate_compare_coefficients(NUM_LINE_POINTS,H2[0])
        self.line7 = visuals.Line(line_data7, parent=self.view_bot4.scene, color=LINE_COLOR_CHOICES[0], width = 2)
        self.line8 = visuals.Line(line_data8, parent=self.view_bot4.scene, color=LINE_COLOR_CHOICES[1], width = 2)
        self.view_bot4.camera = "panzoom"
        self.view_bot4.camera.set_range(x=(0, 1), y=(0, 0.6))         # x and y are used to set the range of the axes

        yaxis4 = scene.AxisWidget(orientation='left', axis_label='Ci', axis_font_size=8, axis_label_margin=30, tick_label_margin=2)
        yaxis4.width_max = 50
        self.grid.add_widget(yaxis4, row=6, col=0)
        yaxis4.link_view(self.view_bot3)
        x_axis5 = scene.AxisWidget(orientation='bottom', axis_label='Time (%)', axis_font_size=8, axis_label_margin=20,tick_label_margin=15)
        x_axis5.stretch = (1,0.2)
        self.grid.add_widget(x_axis5, row=7, col=1)
        x_axis5.link_view(self.view_bot3)
        x_axis6 = scene.AxisWidget(orientation='bottom', axis_label='Time (%)', axis_font_size=8, axis_label_margin=20,tick_label_margin=15)
        x_axis6.stretch = (1,0.2)
        self.grid.add_widget(x_axis6, row=7, col=2)
        x_axis6.link_view(self.view_bot4)

    def updating_canvas(self,x,x2):
        first = generate_compare_coefficients(len(x[0]),x[0])
        second = generate_compare_coefficients(len(x2[0]),x2[0])
        self.line.set_data(first, color=LINE_COLOR_CHOICES[0])
        self.line2.set_data(second, color=LINE_COLOR_CHOICES[1])

        first = generate_compare_coefficients(len(x[1]),x[1])
        second = generate_compare_coefficients(len(x2[1]),x2[1])
        self.line3.set_data(first, color=LINE_COLOR_CHOICES[0])
        self.line4.set_data(second, color=LINE_COLOR_CHOICES[1])

        first = generate_compare_coefficients(len(x[2]),x[2])
        second = generate_compare_coefficients(len(x2[2]),x2[2])
        self.line5.set_data(first, color=LINE_COLOR_CHOICES[0])
        self.line6.set_data(second, color=LINE_COLOR_CHOICES[1])

        first = generate_compare_coefficients(len(x[3]),x[3])
        second = generate_compare_coefficients(len(x2[3]),x2[3])
        self.line7.set_data(first, color=LINE_COLOR_CHOICES[0])
        self.line8.set_data(second, color=LINE_COLOR_CHOICES[1])

        self.canvas.update()

    def optimization_canvas(self, iter, x_range, y_pred, y_std, ei, sample_x, sample_y, new_x, new_y):
        if iter == 3:
            firstpred = generate_optimization(x_range,y_pred)
            self.linepred.set_data(firstpred, color=LINE_COLOR_CHOICES[0], width = 2.2)
        else:
            allpred = generate_optimization(x_range,y_pred)
            self.linepred2.set_data(allpred, color=LINE_COLOR_CHOICES[1], width = 2.2)
        vals = np.column_stack((sample_x, sample_y))
        new_line_ei = generate_optimization(x_range,ei)
        self.line_ei.set_data(new_line_ei, color=LINE_COLOR_CHOICES[0])
        
        self.markers.set_data(vals, edge_color = 'black', face_color = 'black')
        self.themarker.set_data(pos = np.array([[new_x, new_y]]), edge_color = 'green', face_color = 'green')
        # Add the shaded area to the plot CHECK THE LINE BELOW
        vertices = fill_region(x_range,y_pred - 2* y_std,y_pred + 2* y_std)
        self.area._pos = vertices
        self.area._update()

def fill_region(x,y1,y2):
    vertices = np.column_stack((np.concatenate((x, x[::-1])), np.concatenate((y1, y2[::-1]))))
    return vertices
    
def generate_optimization(x_range,values, dtype=np.float32):
    pos = np.empty((len(x_range), 2), dtype=dtype)
    pos[:, 0] = x_range
    pos[:, 1] = values
    return pos

def generate_compare_coefficients(num_points,values, dtype=np.float32):
    x = np.linspace(0,1,num_points)
    pos = np.empty((len(x), 2), dtype=dtype)
    pos[:, 0] = x
    pos[:, 1] = values
         
    return pos