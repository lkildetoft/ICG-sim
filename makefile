default:
  gcc -std=c17 Distrs.c Simulation.c ddHist.c RandChoice.c  main.c -o sim.a -fopenmp O3
