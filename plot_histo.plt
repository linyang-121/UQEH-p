#This script will generate a histogram with gnuplot
#using the fort.32 otuput file.
#The produced plot will be saved as "histo.png" in the directory
#the following command is typed.
#Use the terminal command "gnuplot plot_histo.plt" to 
#generate the desired plot.
     set style fill solid 1.0 border -1
     set title "Distribution of electronic Hamiltonian
     set ylabel "Frequency"
     set xlabel "Value of electronic Hamiltonian"
     set key off
     plot 'fort.32' with boxes
     set terminal pngcairo
     set output "histo.png"
     replot
