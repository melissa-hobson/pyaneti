def create_pars_tr(df_vec,planet_labels=plabels[:nplanets],band=bands[0]):

    # Create parameters vector
    pars_tr = np.zeros(shape=(len(planet_labels), 6))
    rp_vals = []
    for m,p in enumerate(planet_labels):
        pars_tr[m][0] = df_vec[f't0{p}']
        pars_tr[m][1] = df_vec[f'p{p}']
        pars_tr[m][2] = df_vec[f'e{p}']
        pars_tr[m][3] = df_vec[f'omega{p}']
        pars_tr[m][4] = df_vec[f'i{p}']
        pars_tr[m][5] = df_vec[f'arstar{p}']
        if nradius == 1:
            rp_vals.append(df_vec[f'rprstar{p}'])
        else:
            for band in bands:
                rp_vals.append(df_vec[f'rprstar{band}{p}'])

    if nbands == 1:
        us = [df_vec.u1,df_vec.u2]
    else:
        us = [df_vec[f'u1{band}'],df_vec[f'u1{band}']]

    return pars_tr, rp_vals, us

def compute_lc_model_single_band(x,df_vec,planet_labels=plabels[:nplanets],band=bands[0],n_cad=n_cad[0],t_cad=t_cad[0]):

    #Get the parameters
    pars_tr, rp_vals, us = create_pars_tr(df_vec,planet_labels,band=band)
    #
    my_trlab = [0]*len(x)
    f = pti.flux_tr(x, my_trlab, pars_tr.transpose(),rp_vals, us, n_cad, t_cad,1)

    return f

#--------------------------------------------------------------------------------
#                   Light curve plot
#--------------------------------------------------------------------------------

plt.figure(1, figsize=(fsx, fsy/2))
delta_y = df_medians[radius_rstar_name]**2
#extract the data of the current band
for i in range(nbands):
    ix = i == np.array(trlab)
    local_time = lc_time[ix]
    local_flux = lc_flux[ix]
    local_errs = lc_errs[ix]
    if is_jitter_tr:
        local_errs_jit = np.sqrt(lc_errs[ix]**2 + df_medians[f'tr_jitter{bands[i]}']**2)
    else:
        local_errs_jit = local_errs

    # Here I need to create a special trlab in order to separate the different colors
    # Now let us imagine that it works with 1-band
    xmodel = np.linspace(min(local_time)-1, max(local_time)+1,5000)

    lc_model = compute_lc_model_single_band(xmodel,df_medians,band=bands[i],n_cad=n_cad[i],t_cad=t_cad[i])
    lc_residuals = compute_lc_model_single_band(local_time,df_medians,band=bands[i],n_cad=n_cad[i],t_cad=t_cad[i])
    lc_residuals = local_flux - lc_residuals

    plt.plot(local_time,local_flux-delta_y*i,'o',mfc='w',label=bands[i],color=tr_colors[i])
    plt.plot(xmodel,lc_model-delta_y*i,color=tr_colors[i])

    fn = f'{outdir}/{star}_{bands[i]}_lightcurve_model.dat'
    with open(fn,'w') as f:
        f.write('#Time  flux \n')
        for j in range(len(lc_model)):
            f.write(f'{xmodel[j]} {lc_model[j]} \n')

    fn = f'{outdir}/{star}_{bands[i]}_lightcurve_data.dat'
    with open(fn,'w') as f:
        f.write('#Time  flux  eflux eflux_jitter flux_residuals \n')
        for j in range(len(local_time)):
            f.write(f'{local_time[j]} {local_flux[j]} {local_errs[j]} {local_errs_jit[j]} {lc_residuals[j]} \n')

plt.legend()
plt.xlim(lc_time.min()-1,lc_time.max()+1)
plt.xlabel(tr_xlabel)
plt.ylabel('Normalised flux')



fname = outdir+'/'+star+'_lightcurve.pdf'
print(f'Creating  {fname}')
plt.savefig(fname, format='pdf', bbox_inches='tight',dpi=300)
plt.savefig(fname[:-3]+'png', format='png', bbox_inches='tight', dpi=300)
plt.close()

#--------------------------------------------------------------------------------
#                   Phase folded transit plots
#--------------------------------------------------------------------------------

#Create phase-folded plots for all bands, for all planets
for k,p in enumerate(plabels[:nplanets]):

    if not fit_tr[k]: continue

    plt.figure(1, figsize=(fsx, fsy/1.618 + fsy*(nbands-1)*0.3))

    gs = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[fsy/1.618 + fsy*(nbands-1)*0.3, fsy/1.618/2.])
    gs.update(hspace=0.00)
    ax0 = plt.subplot(gs[0],rasterized=is_rasterized) #models and data plots
    ax1 = plt.subplot(gs[1],rasterized=is_rasterized) #residuals

    #Let's find the planet radius that we will use for the next calculations
    radius_pattern = f"rprstar{p}"
    # Extract the first element that contains the string pattern
    radius_rstar_name = next((x for x in list(df) if radius_pattern in x), None)
    #radius_name contains the df name of the first element of radius rp for {planet}
    delta_y = df_medians[radius_rstar_name]**2

    #extract the data of the current band
    for i in range(nbands):

        ix = i == np.array(trlab)
        local_time = lc_time[ix]
        local_flux = lc_flux[ix]
        local_errs = lc_errs[ix]
        if is_jitter_tr:
            local_errs_jit = np.sqrt(lc_errs[ix]**2 + df_medians[f'tr_jitter{bands[i]}']**2)
        else:
            local_errs_jit = local_errs

        if nplanets > 1:
            #compute the models for the planets we are not modeling
            other_planets = [label for label in plabels[:nplanets] if label != p]
            #Now we compute the model for the other planets
            lc_op = compute_lc_model_single_band(local_time,df_medians,other_planets,band=bands[i],n_cad=n_cad[i],t_cad=t_cad[i])
            #Remove it from the data for the given band
            local_flux = local_flux/lc_op

        #Compute the model for the current planet
        xmodel = np.linspace(df_medians[f't0{p}']-1.5*df_medians[f'trt{p}']/24.,df_medians[f't0{p}']+1.5*df_medians[f'trt{p}']/24.,500)
        fmodel = compute_lc_model_single_band(xmodel,df_medians,p,band=bands[i],n_cad=n_cad[i],t_cad=t_cad[i])
        fdata  = compute_lc_model_single_band(local_time,df_medians,p,band=bands[i],n_cad=n_cad[i],t_cad=t_cad[i])
        #other samples
        #compute random samples
        #fsamples = []
        #for j in range(n_samples):
        #    fsamples.append(compute_lc_model_single_band(xmodel,df_samples.iloc[j],p,band=bands[i],n_cad=n_cad[i],t_cad=t_cad[i]))
        xmodel -= df_medians[f't0{p}']
        pmodel = xmodel*24

        #Let us phase fold the light curve
        phase = abs(((local_time-df_medians[f't0{p}'])%df_medians[f'p{p}'])/df_medians[f'p{p}'])
        phase[phase>0.5] -= 1
        phase *= df_medians[f'p{p}']*24.
        indices = abs(phase) <= 1.5*df_medians[f'trt{p}']

        #Let us take care of planets that were not observed with a given band
        if len(local_flux[indices]) == 0: continue

        lc_residuals = local_flux[indices] - fdata[indices]

        ax0.plot(phase[indices],local_flux[indices]-i*1.5*delta_y,mark_tr[i],color='k',alpha=0.05)
        pb, lfb, efb = bin_data(phase[indices],local_flux[indices],local_errs[indices],tbin)
        ax0.errorbar(pb,lfb-i*1.5*delta_y,efb,mfc='w',label=bands[i],color=tr_colors[i],zorder=6,fmt=mark_tr[i])
        ax0.plot(pmodel,fmodel-i*1.5*delta_y,zorder=5,color=tr_colors[i])
        #for j in range(n_samples):
        #    ax0.plot(pmodel,fsamples[j]-i*1.5*delta_y,alpha=0.05,lw=0.5,zorder=4,color=tr_colors[i])

        #Residuals
        pb, lfb, efb = bin_data(phase[indices],lc_residuals,local_errs[indices],tbin)
        ax1.plot(pb,lfb*1e3,marker=mark_tr[i],mfc='w',color=tr_colors[i],zorder=6,alpha=0.75,lw=0)

        #Store the current band for the current planet

        fn = f'{outdir}/{star}{p}_tr_{bands[i]}_model.dat'
        with open(fn,'w') as f:
            f.write('#Time  flux \n')
            for j in range(len(pmodel)):
                f.write(f'{pmodel[j]} {fmodel[j]} \n')

        fn = f'{outdir}/{star}{p}_tr_{bands[i]}_data.dat'
        with open(fn,'w') as f:
            f.write('#Time  flux  eflux eflux_jitter  flux_residuals \n')
            for j in range(len(phase[indices])):
                f.write(f'{phase[indices][j]} {local_flux[indices][j]} {local_errs[indices][j]} {local_errs_jit[j]} {lc_residuals[j]} \n')

    ax0.set_ylim(1-2*delta_y*(nbands),1+delta_y)
    ax0.set_xlim(-1.5*df_medians[f'trt{p}'],1.5*df_medians[f'trt{p}'])
    ax0.legend()
    ax0.set_xticks([])
    ax0.set_ylabel('Flux')
    ax0.ticklabel_format(useOffset=False, axis='y')

    ax1.set_ylabel('Residuals \n (ppt)')
    ax1.set_ylim(-delta_y*1e3,delta_y*1e3)
    ax1.set_xlim(-1.5*df_medians[f'trt{p}'],1.5*df_medians[f'trt{p}'])
    ax1.set_xlabel('Time [hours]')

    fname = f'{outdir}/{star}{p}_tr.pdf'
    print(f'Creating  {fname}')
    plt.savefig(fname, format='pdf', bbox_inches='tight',dpi=300)
    plt.savefig(fname[:-3]+'png', format='png', bbox_inches='tight', dpi=300)
    plt.close()

