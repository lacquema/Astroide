#! /Users/lacquema/Oracle.env/bin/python3
import sys

from PyQt6.QtCore import pyqtSignal
from PyQt6.QtWidgets import QApplication, QSplitter, QMainWindow, QStatusBar

from WidgetParam import WidgetParam
from WidgetPlot import WidgetPlot

class WindowPlot(QMainWindow):

    SignalCloseWindowPlot = pyqtSignal()  # initiation of the closeEvent signal

    def __init__(self, ToolName):
        super().__init__()

        # Window settings
        self.setWindowTitle(ToolName)

        # Splitter widget
        self.Splitter = QSplitter()

        # Parameters widget
        self.WidgetParam = WidgetParam()
        self.Splitter.addWidget(self.WidgetParam)

        # Plotting widgets
        self.WidgetPlots = []

        # Status bar
        self.setStatusBar(QStatusBar(self))

        # Container
        self.setCentralWidget(self.Splitter)


    def add_WidgetPlot(self, plot, events_to_reset_history=None):
        """
        Creates a new WidgetPlot connected to the WidgetParam.
        """
        widget_plot = WidgetPlot(plot)
        self.connect_events_to_reset_history(widget_plot, events_to_reset_history)
        self.WidgetPlots.append(widget_plot)
        self.Splitter.addWidget(widget_plot)
        return widget_plot
    
    def connect_events_to_reset_history(self, widget_plot, events_to_reset_history=None):
        """
        Connects the events to reset the history of the widget_plot.
        """
        if events_to_reset_history is None:
            return
        else:
            for x in events_to_reset_history:
                x.connect(widget_plot.reset_history)
    
    


    # Emission of the CloseEvent signal when the parameter window is closed
    def closeEvent(self, e):
        self.SignalCloseWindowPlot.emit()


if __name__=="__main__":
    app = QApplication(sys.argv) # Application creation
    Window = WindowPlot('test')
    Window.show()
    app.exec() # Application execution