Filter Design and Simulation for Solar Inverters

This repository contains MATLAB scripts for the analysis, design, and simulation of filters used in solar inverter systems. The framework supports evaluating different filter topologies to identify factors that enable the smallest and most efficient designs. The LCL filter is provided as an example case study.

Contents
thesis_filter_sim.m – Main script for filter design and simulation.
compute_thd.m – Function to calculate Total Harmonic Distortion (THD) for assessing power quality.

Capabilities
Analysis of multiple filter topologies for inverter applications.
Evaluation of design trade-offs affecting size and performance.
Example implementation and results for an LCL filter.
THD computation to measure harmonic distortion reduction.

Requirements
MATLAB R2022b or later (earlier versions may also work).
Simulink is not required.

How to Run
Clone or download this repository.
Place thesis_filter_sim.m and compute_thd.m in the same folder.

Open MATLAB.
Run:
thesis_filter_sim


Review the output plots and THD values to compare filter performance.

Purpose
This project is based on my undergraduate thesis, which focused on the design, analysis, and simulation of the smallest possible filters for both single-phase and three-phase solar inverter systems, including hybrid and grid-tied configurations.

Contact
Stalin Learoy Koster
Email: koster263@gmail.com
LinkedIn: linkedin.com/in/stalin-koster-747b61b1
