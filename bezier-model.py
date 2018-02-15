"""
Supplementary material S7 for Modelling Human Hard Palate Variation with Bezier Curves
Rick Janssen, Scott R. Moisik, Dan Dediu

File name: BezierSupp.py
Author: Rick Janssen
E-mail: rick.janssen@mpi.nl
Date last modified: 12-01-2016
Python Version: 2.7.5
Licence: GPL

This code contains two classes: BezierCurve that can be used to model the human hard palate on the
mid-sagittal plane using a Bezier curve, InteractivePlotter that is a helper class that 
BezierCurve uses to draw itself into a matplotlib interactive GUI, and VoidPlotter that can be used
if you just want to use the Bezier model without visualization.

The code as it is presented here fires up a GUI that allows you to visualize the Bezier hard 
palate model. On the bottom we have two radiobuttons and five sliders.

The left radiobuttion selects in how much detail we want to visualize the curve.
    -spline: Only show the Bezier curve itself
    -control: Shows the curve, control points and the (top-level) control vectors
    -bezier: Show all of the above plus the recursion vectors
    
The right radiobutton selects the zoom-level.
    -spline: Zoom in so that the curve is visible.
    -control: Zoom in so that the control vectors are visible as well.
    
More to right, we have five sliders that control the appearance of the curve.
    -smoothness: control the sampling interval the curve is interpolated with.
    -palatal concavity: increases the vertical displacement between velar transition and palatal roof.      
    -alveolar angle: controls the angle of inclination of the alveolar ridge from 180 to 90 degr.
    -alveolar weight: modifies the `magnitude' alveolar angle.
    -palatal fronting shifts the palatal roof more anteriorly for higher values.
    
When the left radiobutton is set to either 'control' or 'bezier', the curve's control points are
visible. We can click and drag the control points the change the appearance of the curve in a more
free-form manner. We can also click on the plot to append more control points to the curve
(inserting is as of yet unsupported). Note that we directly manipulating the control points, the 
parameters sliders at the bottom of the screen no longer function properly. If you want to use
the parameters again, we recommend restarting the program.

This version of the model has been cleaned-up and streamlined for the sake of presentation. For our
experiments, we use a version that has functionality accomodating instantiation in a custom-built
evolutionary algorithm and in a Cython/C++ environment. Interested persons should get in contact
with the author.
"""

from math import copysign
import matplotlib.pyplot as plt
import matplotlib.lines as lines
import matplotlib.cm as cm
from matplotlib.widgets import RadioButtons, Slider

class BezierCurve(object):
    """Contains all the methods to to evaluate a Bezier curve following an implementation of
    De Casteljau's algorithm
    
    Attributes:
        control_points: The curve's control points (i.e., Bernstein coefficients)
        spline_points: The curve's n (x,y)-values. You can change the number of values n by calling
            'update_interval(int)'
        params: Dictionary containing the parameter names and associated values
    """
    
    def __init__(self, plotter):
        """Instantiate a Bezier Curve.
        
        Args:
            plotter: Responsible for plotting graphics to e.g., interactive GUI.
        """        
        self._params = {'fronting': 0.5, 'conc': 0.5, 'angle' : 0.5, 'weigth' : 0.5, 
                       'interval': 100} # default parameters 
        self._control_points = [[-0.2, 0.6], [0.1, 0.6], [0.4,0.75], 
                                [0.55, 0.4125], [0.7, 0.3]] # default control points
        self.plotter = plotter
       
    @property
    def control_points(self):
        return self._control_points
    
    @property
    def spline_points(self):
        return self._spline_points
    
    @property
    def params(self):
        return self._params  
       
    def sample(self, intervals=None, n_x=None, min_x=0, max_x=1, min_y=0, max_y=1):
        """Sample the curve's y-values for specific x-values.
        
        Args:
            intervals: You can give the raw x-values here
            n_x: You can also give the number of equidistant sampling points on the x-axis. Note:
                You have to specify either 'interval' or 'n_x'!
            min_x: The minimum horizontal value of the coordinate space you want to translate to
            max_x: The maximum horizontal value of the coordinate space you want to translate to
            min_y: The minimum vertical value of the coordinate space you want to translate to
            max_y: The maximum vertical value of the coordinate space you want to translate to
            
        Returns:
            intervals: The translated sampling intervals
            y_results: The translated y-values sampled
        """
        
        # Calculate x-intervals if none given
        if intervals == None:
            intervals = [i * ((max_x - min_x) / (float(n_x) - 1)) + min_x for i in xrange(n_x)]
        
        # Calculate the curve's profile
        spline_points = self.calculate()
        
        # nNrmalize sampling points on the x-axis to fall within [0.1]
        (min_x_ext, max_x_ext) = (intervals[0], intervals[-1])        
        (min_x_spline, max_x_spline) = (min(zip(*spline_points)[0]), max(zip(*spline_points)[0]))
        
        sample_intervals = []
        for interval in intervals:
            norm_interval = (interval - min_x_ext) / (max_x_ext - min_x_ext)
            norm_interval = norm_interval * (max_x_spline - min_x_spline) + min_x_spline
            sample_intervals.append(norm_interval)
        
        # Sample and interpolate Bezier curve using normalized x-intervals. The y-values are 
        # between [0,1] 
        paired_points = zip(spline_points[:-1], spline_points[1:])
        
        samples = []
        while len(sample_intervals) != 0:
            try:
                while paired_points[0][1][0] < sample_intervals[0]:
                    paired_points = paired_points[1:]
                else:
                    (x0,y0) = (paired_points[0][0][0], paired_points[0][0][1])
                    (x1,y1) = (paired_points[0][1][0], paired_points[0][1][1])
                                
                    y = y0 + (y1 - y0) * ((sample_intervals[0]  - x0) / (x1 - x0))
                     
                    samples.append([sample_intervals[0] ,y])
                    sample_intervals = sample_intervals[1:]
            except IndexError:
                samples.append(spline_points[-1])
                sample_intervals = sample_intervals[1:] 
        
        # Translate y-values so that they fall within [y_min,y_max] 
        y_samples = zip(*samples)[1]
        (measured_min_y, measured_max_y) = (y_samples[-1], y_samples[0])        
        
        y_results = []
        for i in xrange(len(intervals)):            
            (measured_max_y - measured_min_y)
            normalized = (samples[i][1] - measured_min_y) / (measured_max_y - measured_min_y)
            denormalized = normalized * (max_y - min_y) + min_y
            y_results.append(denormalized)
        
        # return input intervals and translated y-values     
        return (intervals, y_results)       
       
    def calculate(self):
        """Calculate the Bezier Curve's shape. Do this after you have set the curve's parameters
        to see what effect they have.
        """
        
        self._spline_points = [self.control_points[0]] # the curve's raw coordinates
        self.plotter.reset() # init gfx plotter
        
        # Initiate Bezier recursion. Say we sample the curve with 100 sampling points. We iterate 
        # over points in sequential order. For each point, we initiate a recursive procedure where
        # we connect points with vectors (beginning at the top-level control-points. 
        # We recursive draw more vectors between parent-vector on the same point.
        for sample_point in xrange(1, self._params['interval']):
            spline_point = self.descent(self.control_points, sample_point)
            self.spline_points.append(spline_point)
        
        self.spline_points.append(self.control_points[-1])
        
        # draw control points
        for i_line in xrange(len(self.spline_points) - 1):
            self.calc_line(self.spline_points, i_line)                 
        self.plotter.draw_control_points(self.control_points)
                         
        return self.spline_points
    
    def descent(self, vector_points, sample_point, depth=1):                
        """
        Recursively descent a Bezier curve for a given sample-point.
        
        Args:
            vector_points: A list of 2D-points defining the vectors a particular depth of the curve.
            interval_point: A float defining which sampling point we are descending.
            depth: Keeps track of the depth we're at.
            
        Returns:
            vector_points: A list of 2D-points that define the Bezier curve.
        """     
        
        # Use recursion as long as we can construct pairs of points between vectors.
        if len(vector_points) > 1:
            
            # First some variables for drawing the curve
            c_code = ((float(depth) - 1) / (len(self.control_points) - 2))
            zorder = depth - len(self.control_points) - 1
        
            if __name__ == '__main__':
                color = cm.cool(c_code) # @UndefinedVariable
            else:
                color = None
            
            spline_draw = depth==1 and not self.plotter.draw_level=='spline' and sample_point==1
            draw = True if self.plotter.draw_level=='bezier' or spline_draw else False                        
            # Then the actual calculations
            recursion_points = []
            
            for i_line in xrange(len(vector_points) - 1):
                (x1,y1,x2,y2) = self.calc_line(vector_points, i_line, draw, color, zorder)
                
                dx = (x2 - x1) / self._params['interval']
                dy = (y2 - y1) / self._params['interval']       
                x = x1 + dx * sample_point
                y = y1 + dy * sample_point           
                recursion_points.append([x,y])
                
            # Do the same stuff with the newly calculated points and recursively return appended.
            return self.descent(recursion_points, sample_point, depth+1)
        # If we ran out of vectors, simply return the endpoint of the curve.
        else:
            return vector_points[0]

    def calc_line(self, points, i, draw=True, color='mediumblue', zorder=0):
        """Fetch the ith vector from a list of 2D-points.
        
        Args:
            points: A list of 2D-points.
            i: Which vector to fetch.
            draw: Do we want to visualize the vector?
        
        Returns:
            vector: A tuple (x0,y0,x1,y1)
        
        """
        vector = (points[i][0], points[i][1], points[i+1][0], points[i+1][1])

        if draw:
            self.plotter.draw_line_segments(vector, zorder, color)
            
        return vector
    
    def update_params(self):
        """Update control points. Call this after having set the parameters.
        """        
        
        a_param = self._params['angle']
        c_param = self._params['conc']
        f_param = self._params['fronting']
        w_param = self._params['weigth']
        
        fronting = self.control_points[1][0] + (self.control_points[-1][0] - self.control_points[1][0]) * f_param
        self.control_points[2] = [fronting, self.control_points[1][1] + (self.control_points[0][1] - self.control_points[-1][1]) * c_param]
        
        palate_x = (self.control_points[-1][0] - self.control_points[1][0]) * (1 - a_param)
        palate_x = self.control_points[-1][0] - palate_x * w_param
        
        palate_y = (self.control_points[-1][1] - self.control_points[2][1]) * a_param
        palate_y = self.control_points[-1][1] - palate_y * w_param
    
        self.control_points[3] = [palate_x, palate_y]   
    
    
    #Below are the parameter update functions (hooks for e.g., GUIs, optimisation algorithms, etc.)
    def update(self, update_params=True, update_curve=True):
        """Updates the Bezier curve
        
        Args:
            update_params: Do you want to update the parameters? (set to 'False' for manual 
                adjustment).
            update_curve: Do you want to update the curve profile? (set to 'False' to adjust
                multiple parameters before finalizing)
        """
        
        if update_params:
            self.update_params()
            
        self.sample(n_x=100)
    
    def update_fronting(self, value):
        """ Update the fronting parameter and update curve.
        
        Args:
            value: Set fronting to this value
        """
        
        self._params['fronting'] = float(value)
        self.update()
                 
    def update_conc(self, value):
        """ Update the concavity parameter and update curve.
        
        Args:
            value: Set concavity to this value
        """     
           
        self._params['conc'] = float(value)
        self.update()
        
    def update_angle(self, value):
        """ Update the angle parameter and update curve.
        
        Args:
            value: Set angle to this value
        """  
              
        self._params['angle'] = float(value)
        self.update()
        
    def update_weigth(self, value):
        """ Update the fronting parameter and update curve.
        
        Args:
            value: Set weight to this value
        """   
             
        self._params['weigth'] = float(value)
        self.update()
            
    def update_interval(self, value):
        """Set the sampling interval and update curve.
        
        Args:
            value: Set sampling interval to this value
        """   
             
        value = int(value)
        
        try:
            if value % int(value) != 0:
                self.intervalSlider.set_val(int(value))
            else:
                self._params['interval'] = value
                
            self.update()
        except ZeroDivisionError:
            pass

class InteractivePlotter():
    """Used to plot Bezier curve to matplotlib window.
    
    Properties:
        curve: The Bezier curve instance which to plot
        draw_level: How deep top we want to plot (level 1 for only control vectors)
    """
    
    def __init__(self, ax):
        """Instantiate plotter
        
        Args:
            ax: A matplotlib axes instance to which to plot.
        """
        
        self.ax = ax
        (self.boundaries,self.increments) = ([],[])
        self._draw_level = 'control'
        self.zoom = self.draw_level
    
    @property
    def curve(self):
        return self._curve
    
    @curve.setter
    def curve(self, curve):
        self._curve = curve
        
    @property
    def draw_level(self):
        return self._draw_level
    
    def draw_line_segments(self, vector, zorder, color):
        """Draws line segments (e.g., control vectors and curve segments.
        
        Args:
            vector: A tuple (x0,y0,x1,y1)
            zorder: z-position of the line segment.
            color: Color of the line segment
        """
        
        (x1,y1,x2,y2) = vector
        lw = 4 if zorder == 0 else 1.5
        line1 = lines.Line2D([x1, x2],[y1, y2], color=color, marker='', ms=8, zorder=zorder,
                                           mec='None', lw=lw)
        self.ax.add_line(line1)
        
    def draw_control_points(self, points):
        """Draws Bezier curve control points
        
        Args:
            points: A tuple (x0,y0,x1,y1)
        """
                
        if self.draw_level != 'spline':
            zpoints = zip(*points)
            self.ax.plot(zpoints[0], zpoints[1], marker='o', ls='o', ms=20, mfc='mediumblue', 
                    mec='None', picker=10)
            
        self.zoom_and_draw()     
    
    def reset(self):
        """ Resets the axes to initial configuration.
        """
        
        self.ax.cla()
        self.ax.minorticks_on()
        self.ax.grid(which='both')
        
    def zoom_and_draw(self):
        """Set the appropriate zoom level.
        """
    
        if self.zoom == 'spline':
            boundaries = self.get_axis_increments(self.curve.spline_points)
        elif self.zoom == 'control':
            boundaries = self.get_axis_increments(self.curve.control_points)
           
        self.ax.axis(boundaries)
            
        plt.draw() 
        
    def get_axis_increments(self, points):
        """Calculate axis limits based on points displayed.
        
        Args:
            points: A list of 2D-points
        """
                
        points = (zip(*points)[0], zip(*points)[1])
        extremes = [(min(axis), max(axis)) for axis in points]
        increments = [(axis[1] - axis[0]) * 0.05 for axis in extremes]
        
        boundaries = []
        for i in xrange(2):
            for j in xrange(2):            
                boundary = extremes[i][j] + copysign(increments[i], j - 1)
                boundaries.append(boundary)
                  
        return boundaries  
        
    def update_level(self, label):
        """This sets the recursion depth to plot.
        
        Args: 
            label: 'control', 'bezier' or 'spline'
        """
        
        self.draw_level = label
        self.curve.update()
        
    def update_zoom(self, label):
        """This sets the zoom level.
        
        Args: 
            label: 'control' or 'spline'
        """
                
        self.zoom = label
        self.zoom_and_draw()
        
class voidPlotter():
    """Use this plotter is you don't want to visualize anything but just want to use to Bezier
    model without graphical slowdown.
    """
    
    def draw_line_segments(self):
        pass
    
    def reset(self):
        pass
    def draw_control_points(self):
        pass
    
    def draw_level(self):
        pass
        
#Instantiate and setup up matplotlib GUI
if __name__== "__main__":
    # instantiate matplotlib
    fig = plt.figure(figsize=(18, 12))
    min_canvas_y = 0.15
    ax = fig.add_axes([0.03, min_canvas_y, 0.95, 0.82])    
    plotter = InteractivePlotter(ax)
    
    # instantiate Bezier curve
    bezier = BezierCurve(plotter)
    plotter.curve = bezier
    bezier.sample(min_y=0, max_y=1, n_x=10)
    
    #instantiate GUI sliders
    axFronting = fig.add_axes([0.22, 0.01, 0.7, 0.01])
    frontingSlider = Slider(axFronting, 'fronting', 0, 1, bezier.params['fronting'])
    frontingSlider.on_changed(bezier.update_fronting)
    
    axPalateConc = fig.add_axes([0.22, 0.025, 0.7, 0.01])
    concSlider = Slider(axPalateConc, 'palate concavity', 0, 1, bezier.params['conc'])
    concSlider.on_changed(bezier.update_conc)
    
    axAngle = fig.add_axes([0.22, 0.04, 0.7, 0.01])
    angleSlider = Slider(axAngle, 'alveolar angle', 0, 1, bezier.params['angle'])
    angleSlider.on_changed(bezier.update_angle)
    
    axWeigth = fig.add_axes([0.22, 0.055, 0.7, 0.01])
    weightSlider = Slider(axWeigth, 'palatal weigth', 0, 1, bezier.params['weigth'])
    weightSlider.on_changed(bezier.update_weigth)
    
    axInterval = fig.add_axes([0.22, 0.1, 0.7, 0.01])
    intervalSlider = Slider(axInterval, 'smoothness', 0, 100, valinit=bezier.params['interval'], color='r')
    intervalSlider.on_changed(bezier.update_interval)
    
    # instantiate radio buttons
    axShow = fig.add_axes([0.02, 0.02, 0.05, 0.09])
    radio_options = ('control', 'bezier', 'spline')
    showRadio = RadioButtons(axShow, radio_options, active=radio_options.index(plotter.draw_level))
    showRadio.on_clicked(plotter.update_level)
    
    axZoom = fig.add_axes([0.075, 0.02, 0.05, 0.09])
    zoomOptions = ('control', 'spline')
    zoomRadio = RadioButtons(axZoom, zoomOptions, zoomOptions.index(plotter.draw_level))
    zoomRadio.on_clicked(plotter.update_zoom)  
    
    # Below are the control point handlers
    picked = None # keeps track if we picked up a control point 
    
    def onpick(pickEvent):
        global picked
        picked = pickEvent.ind
    
    def onpress(mouseEvent):
        global picked
        if picked == None:
            point = [mouseEvent.xdata, mouseEvent.ydata]
            max_y = fig.canvas.get_width_height()[1] * min_canvas_y           
            if point != [None, None] and mouseEvent.y > max_y:
                bezier.control_points.append(point)
    
    def onrelease(mouseEvent):
        global picked
        if picked != None:
            bezier.control_points[picked] = (mouseEvent.xdata, mouseEvent.ydata)         
        picked = None
        
        point = [mouseEvent.xdata, mouseEvent.ydata]
        max_y = fig.canvas.get_width_height()[1] * min_canvas_y

        if point != [None, None] and mouseEvent.y > max_y:
            bezier.update(update_params=False)

    # connect buttons to fucntions
    fig.canvas.mpl_connect('button_release_event', onrelease)
    fig.canvas.mpl_connect('button_press_event', onpress)
    fig.canvas.mpl_connect('pick_event', onpick)      
        
    # finally, show the GUI    
    plt.show()