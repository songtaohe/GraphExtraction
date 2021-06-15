# Usage
```bash
usage: main.py [-h] -input INPUT -output OUTPUT [-threshold THRESHOLD]
               [-gaussian_sigma GAUSSIAN_FILTER_SIGMA]
               [-grey_closing GREY_CLOSING]
               [-connected_component_min_size CONNECTED_COMPONENT_MIN_SIZE]
               [-spurs_min_size SPURS_MIN_SIZE]
               [-connect_deadends_threshold CONNECT_DEADENDS_THRESHOLD]
               [-polyline_simplify_threshold POLYLINE_SIMPLIFY_THRESHOLD]

optional arguments:
  -h, --help            show this help message and exit
  -input INPUT          Input segmentation file
  -output OUTPUT        Output file prefix
  -threshold THRESHOLD  Threshold for binarization (0-255)
  -gaussian_sigma GAUSSIAN_FILTER_SIGMA
                        The sigma (pixels) of the gaussian filter. Set this to
                        0 to disable gaussian filter.
  -grey_closing GREY_CLOSING
                        The window size (pixels) for grey_closing
  -connected_component_min_size CONNECTED_COMPONENT_MIN_SIZE
                        Remove connected components whose total lengths are
                        less than [connected_component_min_size] pixels
  -spurs_min_size SPURS_MIN_SIZE
                        Remove short edges (deadends) whose lengths are less
                        than [spurs_min_size] pixels
  -connect_deadends_threshold CONNECT_DEADENDS_THRESHOLD
                        Deadends (terminal vertices) within
                        [connect_deadends_threshold] pixels will be connected
  -polyline_simplify_threshold POLYLINE_SIMPLIFY_THRESHOLD
                        We use Douglas Peucker algorithm to simplify the
                        polylines. A larger threshold gives more aggressive
                        simplification
```

# Example
```bash
python main.py -input examples/0.png -output test
```