#
# Julia script to calculate SiRi calibration coefficients
# It takes as input the Mo96_siri_peaks.csv file,
# which is output from the "clicking" script peaks2D.C
# Version for 96Mo, exp. 96Mo(p,p') at OCL in March 2019
# Cecilie, 28 March 2019
#
using FileIO
using CSVFiles, DataFrames

df = DataFrame(load("Mo96_siri_peaks.csv"))

for row in CSV.File(file)
    println("a=$(row.a), b=$(row.b), c=$(row.c)")
end

