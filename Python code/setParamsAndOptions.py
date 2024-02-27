import numpy as np

# Create a dictionary for options
Options = {}

# Period to analyze
Options["PeriodToAnalyze"] = [0, np.inf]

# Save MUA option
Options["SaveMUA"] = 0

# LogMUA options
Options["LogMUA"] = {
   "FreqBand": [200, 1500],
   "SmoothingWindow": 0.040
}

# LFP options
Options["LFP"] = {
   "SmoothingWindow": Options["LogMUA"]["SmoothingWindow"],  # Use LogMUA's smoothing window
   "SubsamplingRate": Options["LogMUA"]["FreqBand"][0],
   "HPFCutOffFreq": 0.2,
   "xWavScale": 1e6
}

# UpTimeLags options
Options["UpTimeLags"] = {
   "MaxAbsTimeLag": 0.250
}

# LFPMUAplot options
Options["LFPMUAplot"] = {
   "MUArange": [-0.75, 3.5],
   "LFPrange": [-400, 400]
}

# MUARef options
Options["MUARef"] = {
   "Method": "LeftTail",
   "MinDensity": 0.05,
   "MUAShift": 0.35
}

# PrintDev option
PrintDev = '-dpdf'

# Define colors and palettes
Red = np.array([1, 0, 0])
Green = np.array([0, 1, 0])
Blue = np.array([0, 0, 1])
Yellow = Red + Green
YellowOrange = 0.75 * Yellow + 0.25 * Red
Orange = 0.5 * Yellow + 0.5 * Red
BluePurple = 0.75 * Blue + 0.25 * Red

# Define a custom function for graded colormap (assuming you'll define it later)
RedBlueCB = gradedColormap(Blue, Red)  # Replace with implementation if needed
