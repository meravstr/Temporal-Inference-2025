# Inferring population spiking rate from wide-field calcium imaging
Tools for inferring population spiking rate from wide-field (mesoscale) fluorescence. 
The implementation here is in MATLAB; a Python implementation can be provided/posted quickly upon request.
The deconvolution implemented here assumes each recorded fluorescence trace aggregates fluorescence from many active neurons, as in wide-field imaging, and hence reflects their summed activity. Therefore, the deconvolution retrieves the population rate per pixel or trial, not individual spikes.
The repository has one main inference script that utilizes one chosen deconvolution method out of five deconvolution methods provided, plus supporting tools that include:
•	A cross-validated search for each method’s free parameter.
•	Splitting prior-inference long fluorescence traces for fast deconvolution.
•	Concatenating post-inference long spiking rate back to original recording times.
•	A helper for finding your indicator’s calcium decay.
Every line you need to change for your own data is marked % <<< UPDATE HERE.
These tools accompany our paper on mesoscale activity inference (see https://www.biorxiv.org/content/10.1101/2025.08.13.670164v2). 
However, this repository is self-guided, so you don't need to read the paper to use it. The paper can give you a better idea of which method to use in which scenario and what defines “good inference”.
________________________________________
# Quick start
1.	Clone or download the repository and add it to your MATLAB path.
2.	Open run_inference.m and run it as is. It loads the included example (Musall_G6sData.mat, GCaMP6s, dorsal cortex, 30 Hz) and plots the fluorescence next to the inferred rate.
3.	To use your own data, edit the HEADER section of run_inference.m:
4.	
o	Load your data as a time × pixels/trials matrix of ΔF/F named dffed_fluor.
o	Set recording_rate (Hz).
o	Set gamma, the calcium decay per time bin (see calcium_decay_finder.m for help or write to us; your gamma should typically be a number around 0.8- 0.98)
o	Choose method (1–5). We strongly recommend 1 or 4. If you wish to skip parameter search, use 5.
o	Set param_known = true and give param_value if you know the right parameter. Otherwise, leave it false and the script will search for one.
That’s it. The script handles splitting (”chunking”) the parameter search, inference, concatenating, and plotting.
________________________________________
# What’s in the repository
File	What it does
run_inference.m	Start here. A step-by-step script covering the whole pipeline.
calcium_decay_finder.m	Finds gamma for your indicator at your recording rate.
run_method.m	Runs any of the five methods in one call: [r, r_idx_start, r0, beta0] = run_method(y, gamma, method, param).
search_best_param_oddeven.m	Chooses a method’s free parameter by odd/even cross-validation.
chunk_trace.m	Splits long traces into overlapping chunks for fast inference.
stitch_chunks.m	Joins the inferred chunks back into full traces, aligning their baselines.
rebuild_calcium.m	Rebuilds the inferred calcium from an inferred rate, for comparison with fluorescence.
convar.m, dynbin_wstop.m, firdif.m, fft_wiener.m, lucric.m	The five inference methods.
Musall_G6sData.mat	Example data: GCaMP6s, dorsal cortex, 30 Hz, many trials.
ClancyFluorSpikes.mat	Example data: Gcamps6f, v1, 20hz, 5 20min trials
________________________________________
# The pipeline
run_inference.m runs these steps in order:
1.	Header. Load your data and set gamma, the method, and its parameter or search range.
2.	Chunking (chunk_trace.m). Deconvolution time grows quickly with trace length, so traces longer than chunk_size (default 600 samples) are split into chunks that overlap by a quarter of their length. Short, trial-structured data passes through unchanged.
3.	Parameter search (search_best_param_oddeven.m). Skipped if you already know the parameter. Otherwise, see below.
4.	Inference (run_method.m, written out in full in the script so you can see how each method is called).
5.	Stitching (stitch_chunks.m). Joins chunked results back into full traces. Since fluorescence doesn’t include information about baseline firing rates, each chunk is shifted to match global baseline activity, ensuring positive population spiking rates.
6.	Plots. All traces with their average, plus one example trial or pixel. When the data was not split for inference, the example also shows the rebuilt calcium over the fluorescence.
________________________________________
Choosing a method
Method	                          Parameter	  Meaning of the parameter	                  Default search range	    Inferred spiking rate from t = (r_idx_start)
1	Continuously-Varying (Convar),  lambda	    Smoothness penalty; larger is smoother	    logspace(-4, 3, 50)	      (2)
2	Dynamically-Binning	(Dynbin)    lambda	    Smoothness penalty; larger is smoother	    logspace(-4, 3, 50)	      (2)
3	First-Differences (Firdif)	    smt	        Smoothing window (samples)	                1-50	                    (2)
4	Wiener-Filter (Wiener/fft)	    k	          Related to inverse SNR; larger is smoother	logspace(-3, 4, 50)	      (1)
5	Lucy-Richardson (Lucy/Lucric)	  iter	      Number of iterations; 	                    1-29 (default value: 10)	(1)
________________________________________
Finding gamma
gamma is the fraction of calcium signal left after one time bin. It depends on both the indicator and your recording rate. Open calcium_decay_finder.m, set your recording rate and indicator, and run it:
Indicator	Reference gamma	At
GCaMP6f	0.97	40 Hz
GCaMP6s	0.95	10 Hz 
For any other indicator, choose 'custom' and enter a gamma you know along with the rate it was measured at. The conversion between rates is
gamma_new = gamma_ref ^ (rate_ref / rate_new)
________________________________________
Choosing the parameter
When param_known = false, search_best_param_oddeven.m follows Jewell & Witten (2018):
1.	Each trace is split into odd and even time points.
2.	One half is deconvolved and turned back into calcium, which is used to predict the other half’s fluorescence.
3.	This is repeated for every candidate parameter.
4.	The chosen value is the smoothest one whose error stays within one standard deviation of the minimum error.
A plot of the error curve shows the choice. By default, the search uses up to 100 randomly chosen traces (search_maxp in run_inference.m).
If you process many recordings from the same preparation, run the search once, then set param_known = true with that value for the rest.
________________________________________
# Data format
•	Input: a T × N matrix of ΔF/F, where rows are time points, and columns are pixels or trials. Each column is inferred independently. The script rescales the data before inference, since fluorescence units are arbitrary.
•	Output: r_full, a (T − r_idx_start + 1) × N matrix of the inferred population rate, in the same column order as your input.
________________________________________
# Citation
If you use this repository, please cite:
Temporal Deconvolution of Mesoscale Recordings
If you used Continuously-Varying as your inference method of choice, also cite:
An Analytic Framework for Inferring Population Dynamics from Aggregated Calcium Fluorescence
The example data comes from:
Musall, S., Kaufman, M. T., Juavinett, A. L., Gluf, S., & Churchland, A. K. (2019). Single-trial neural dynamics are dominated by richly varied movements. Nature Neuroscience, 22, 1677–1686.
Clancy, K. B., Orsolic, I. & Mrsic-Flogel, T. D. (2019), ‘Locomotion-dependent remapping of
distributed cortical networks’, Nature Neuroscience
________________________________________
Contact
mstern@rockefeller.edu
This helper
Was partly automatically produced from the scripts
