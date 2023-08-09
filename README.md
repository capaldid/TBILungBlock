# TBILungBlock
Tungsten filled 3D printed lung blocks for total body irradiation

This repository is an implementation of "[Tungsten filled 3D printed lung blocks for total body irradiation.](https://capaldi.ucsf.edu)" submitted to the _Physics in Medicine and Biology_ by Dante Capaldi and colleagues.

<p align="center">
  <img width="720" height="348" src="https://github.com/capaldid/TBILungBlock/blob/main/3DTBILungBlocks.jpg">
</p>

## Overview of Matlab Script

3D printed lung blocks, filled with tungsten BBs, are created via a Matlab script, where the input is a Radiation Therapy (RT) Plan DICOM file.  The output is STL files used for 3D printing.

<p align="center">
  <img width="720" height="348" src="https://github.com/capaldid/TBILungBlock/blob/main/MatlabScriptOverview.jpg">
</p>

## Requirements

- MATLAB 2022b or newer
- Radiation Therapy (RT) plan DICOM file (with lung blocks drawn)


## Files and Documentation

The code is provided for lung blocks creation.

```
src\LungBlockGeneration_final.m
```

## Overview of Python Script

3D printed lung blocks for extended field standing TBI procedure, courtesy of J. Schulz.

## Requirements

- Python 3.9

## Files and Documentation

The code is provided for standing lung blocks creation.

```
src\LungBlockGeneration_StandingPython.py
```

## License

MIT License

Copyright 2023 Dante Capaldi

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
