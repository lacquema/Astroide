#! ../bin/python3

# Import necessary packages
import sys
import os
from PyQt6.QtWidgets import QMainWindow, QVBoxLayout, QLabel, QStatusBar, QWidget, QApplication, QPushButton
from PyQt6.QtCore import Qt 

class WindowMenuClass(QMainWindow):
    def __init__(self):
        super().__init__()

        # Get the directory path of the current file
        self.DirPath = os.path.dirname(__file__)

        # Configure the main window
        self.setWindowTitle("Astroide 1.0")
        self.setFixedSize(850, 450)
        self.setStyleSheet(f"background-image: url({self.DirPath}/Items/LoadingBackground.png)") 

        # Initialize the main layout
        Layout = QVBoxLayout()
        Layout.addSpacing(15)

        # Add "New Simulation" button
        self.BtnNew = self.create_button('New Simulation', Layout)
        self.BtnNew.setEnabled(False)
        self.BtnNew.setStyleSheet('QPushButton {background-color: grey; color: grey; font: italic 15px;}')

        # Add "Continue" button
        Layout.addSpacing(20)
        self.BtnContinue = self.create_button('Continue', Layout)
        self.BtnContinue.setEnabled(False)
        self.BtnContinue.setStyleSheet('QPushButton {background-color: grey; color: grey; font: italic 15px;}')

        # Add "Analyse" button
        Layout.addSpacing(90)
        self.BtnAnalyse = self.create_button('Analyse', Layout)

        # Add credits in the status bar
        self.add_status_bar()

        # Set the central widget with the layout
        Container = QWidget()
        Container.setLayout(Layout)
        self.setCentralWidget(Container)

        # Display the window
        self.show()

    def create_button(self, text, layout):
        """
        Helper method to create a styled button and add it to the layout.
        """
        button = QPushButton(text)
        button.setFixedSize(200, 40)
        button.setStyleSheet('QPushButton {background-color: grey; color: white; font: italic 15px;}')
        layout.addWidget(button, alignment=Qt.AlignmentFlag.AlignHCenter)
        return button

    def add_status_bar(self):
        """
        Adds a status bar with credits to the main window.
        """
        StatusBar = QStatusBar(self)
        StatusBar.addWidget(QLabel(' Version 1.0, 2025, IPAG, Hervé Beust, Antoine Lacquement'))
        self.setStatusBar(StatusBar)

if __name__ == "__main__":
    # Create and execute the application
    app = QApplication(sys.argv)
    LoadWin = WindowMenuClass()
    sys.exit(app.exec())