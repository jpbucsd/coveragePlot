#CoveragePlot.py

CoveragePlot was produced by Jordan Burkhardt to reproduce a coverage plot from the paper Comprehensive mapping of the Helicobacter pylori NikR regulon provides new insights in bacterial nickel responses by Vannini, et al. for a school project.
Coverage plot takes a tag directory, complete with bedgraph files produced by homer, as well as the starting and ending location in a bacterial genome and produces coverage plots for that area.

Usage:
python coveragePlot.py -t <tag directory> -d <optional files where each line contains a tag directory name, and the name of the dataset separated by tab> -o <output directory> -s <starting location in genome> -e <ending location in genome> -l <optional flag to add start and stop locations to filenames>
