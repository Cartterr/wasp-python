	real dtdp(61),dtdpsh(61),dtdp_df(61)
	real dist(61),dist_pkp(61)
	data dtdp/8.89,8.84,8.79,8.75,8.69,8.63,8.57,8.51,8.44,
     *            8.37,8.30,8.24,8.17,8.10,8.03,7.96,7.88,7.81,
     *            7.73,7.65,7.58,7.51,7.44,7.37,7.30,7.22,7.15,
     *            7.08,7.00,6.93,6.86,6.80,6.73,6.67,6.60,6.53,
     *            6.47,6.40,6.33,6.26,6.19,6.10,6.01,5.93,5.86,
     *            5.79,5.72,5.64,5.56,5.48,5.40,5.33,5.25,5.16,
     *            5.07,4.98,4.90,4.84,4.78,4.73,4.70/
	data dtdpsh/15.80,15.70,15.60,15.55,15.50,15.45,15.35,
     *              15.25,15.15,15.00,14.90,14.80,14.70,14.55,
     *              14.45,14.30,14.25,14.15,14.00,13.90,13.80,
     *              13.70,13.60,13.50,13.40,13.30,13.20,13.15,
     *              13.00,12.90,12.80,12.70,12.60,12.50,12.35,
     *              12.25,12.10,12.00,11.90,11.80,11.70,11.50,
     *              11.40,11.25,11.15,11.00,10.90,10.80,10.70,
     *              10.50,10.40,10.30,10.15,10.00, 9.90, 9.70,
     *              9.60,9.40,9.30,9.20,8.90/
       data dtdp_df/1.94,1.93,1.93,1.93,1.93,1.93,1.93,1.92,1.92,
     *              1.92,1.92,1.91,1.91,1.90,1.90,1.89,1.88,1.87,
     *              1.86,1.85,1.84,1.83,1.81,1.80,1.78,1.76,1.74,
     *              1.71,1.69,1.66,1.62,1.59,1.55,1.52,1.47,1.43,
     *              1.39,1.34,1.29,1.24,1.19,1.13,1.08,1.02,0.97,
     *              0.91,0.85,0.79,0.73,0.68,0.62,0.55,0.49,0.43,
     *              0.37,0.31,0.25,0.19,0.12,0.06,0.00/
	data dist/30.0,31.0,32.0,33.0,34.0,35.0,36.0,37.0,38.0,
     *            39.0,40.0,41.0,42.0,43.0,44.0,45.0,46.0,47.0,
     *            48.0,49.0,50.0,51.0,52.0,53.0,54.0,55.0,56.0,
     *            57.0,58.0,59.0,60.0,61.0,62.0,63.0,64.0,65.0,
     *            66.0,67.0,68.0,69.0,70.0,71.0,72.0,73.0,74.0,
     *            75.0,76.0,77.0,78.0,79.0,80.0,81.0,82.0,83.0,
     *            84.0,85.0,86.0,87.0,88.0,89.0,90.0/
	data dist_pkp/120.0,121.0,122.0,123.0,124.0,125.0,126.0,
     *          127.0,128.0,129.0,130.0,131.0,132.0,133.0,134.0,
     *          135.0,136.0,137.0,138.0,139.0,140.0,141.0,142.0,
     *          143.0,144.0,145.0,146.0,147.0,148.0,149.0,150.0,
     *          151.0,152.0,153.0,154.0,155.0,156.0,157.0,158.0,
     *          159.0,160.0,161.0,162.0,163.0,164.0,165.0,166.0,
     *          167.0,168.0,169.0,170.0,171.0,172.0,173.0,174.0,
     *          175.0,176.0,177.0,178.0,179.0,180.0/
