# The-Network-Effects-of-Fiscal-Adjustments
Matlab code to replicate the baseline results of the paper. 

The repository contains:
<ol type="i">
  <li><i>Replication_Average_Effects.m:</i> Code to replicate all the results of the paper. 
   <ol>
     <li>The Baysian MCMC code is fully vectorized, however, it still takes quite a few minutes to be run locally (I have a MacBook Pro 15 inches - Processor: 2.3 GHz 8-Core Intel Core i9 - RAM: 16 GB 2400 MHz DDR4.</li>
     <li>The Placebo Simulation repeats 500 baseline Bayesian MCMC simulations. The code is parallelized to speed up performance. My computer took 1h30 minutes to run the code. </li>
    </ol>
  </li>
  <li><i>Functions</i>: This directory contains all Matlab functions to replicate the results. All functions are called by the main code 
      "Replication_Average_Effects.m".</li>
  <li><i>Data</i>: Contains all data required to replicate the results. </li>
</ol>
