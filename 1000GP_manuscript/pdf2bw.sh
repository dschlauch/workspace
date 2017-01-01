#!/bin/bash
gs -sOutputFile=bw_fig.pdf -sDEVICE=pdfwrite -sColorConversionStrategy=Gray -dProcessColorModel=/DeviceGray -dCompatibiltyLevel=1.4 -dAutoRotatePages=/None -dNOPAUSE -dBATCH <(pdftops -level3sep $@ -)
