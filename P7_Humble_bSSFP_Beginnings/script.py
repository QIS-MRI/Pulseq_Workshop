import numpy as np
import pypulseq as pp

# -----------------------
# System & sequence setup
# -----------------------

# Create system limits / hardware parameters
system = pp.Opts(B0=2.89, 
                 max_grad=30, # max gradient amplitude -> confirm with tutor? 
                 grad_unit="mT/m",  
                 max_slew=100, 
                 slew_unit="T/m/s",
                 rf_ringdown_time=30e-6, # time from peak to end of RF (and not block)
                 rf_dead_time=100e-6, # time needed from end to RF1 to start of RF2 
                 adc_dead_time=10e-6, # time needed from end to Acq1 to start of Acq2  
                 adc_raster_time=100e-9, # sampling time vs dwell time? -> basic time unit (typical time in nanos)
                 rf_raster_time=1e-6, # basic rf time unit (typical time in micros)  
                 grad_raster_time=10e-6, # basic gradient time unit
                 block_duration_raster=10e-6 # basic gradient time unit
                 )

seq = pp.Sequence(system=system)

# -----------------------
# Sequence parameters
# -----------------------


name = 'humble_bssfp'

fov = 256e-3  # Define FOV and resolution [m]
delta_k = 1/fov # [1/m]
Nx = 256
Ny = 256
flip_angle_alpha = 30 * 2 * np.pi / 360
# N_dummy_rep = 1 # number of k-space lines (will be computed 
# adc_dur = 500e-6   # seconds

adc_dur = 880e-6
rf_dur = 1000e-6    # seconds
# TR = 5000e-3       # seconds
TE = 30e-3         # seconds
slice_thickness = 5e-3
time_bw = 2



# -----------------------
# Create RF pulse
# -----------------------


rf_excitation, gz_slicesel_pos, _ = pp.make_sinc_pulse(
    flip_angle=flip_angle_alpha/2,
    duration=rf_dur,
    slice_thickness=slice_thickness,
    system=system,
    use='excitation',
    return_gz=True,
    time_bw_product=time_bw
    )


rf_pos, _, _ = pp.make_sinc_pulse(
    flip_angle=flip_angle_alpha,
    duration=rf_dur,
    slice_thickness=slice_thickness,
    system=system,
    use='excitation',
    return_gz=True,
    time_bw_product=time_bw
    )

rf_neg, _, _ = pp.make_sinc_pulse(
    flip_angle=-flip_angle_alpha,
    duration=rf_dur,
    slice_thickness=slice_thickness,
    system=system,
    use='excitation',
    return_gz=True,
    time_bw_product=time_bw
    )

gz_neg = pp.make_trapezoid(
    channel='z',
    system=system,
    area=-gz_slicesel_pos.area/2,
)


gx_readout = pp.make_trapezoid(channel='x', 
                               flat_area=Nx * delta_k, 
                               # flat_time=3.2e-3, 
                               flat_time=adc_dur, 
                               system=system)

gx_neg = pp.make_trapezoid(channel='x', 
                               area=-gx_readout.area / 2, 
                               system=system)

# -----------------------
# ADC event (readout)
# -----------------------

# Calculate delay before ADC:
# delay = TE - RF duration / 2 - RF ringdown time
adc_delay = TE - (rf_dur / 2) - system.rf_ringdown_time

# Safety check
if adc_delay < 0:
    raise ValueError(f"ADC delay is negative: {adc_delay:.6f} s. Check TE and RF settings.")

adc = pp.make_adc(
    num_samples=Nx,
    duration=adc_dur,
    delay=adc_delay,
    system=system
)

# -----------------------
# Delay to complete TR
# -----------------------

# delayTR = TR - pp.calc_duration(rf_pos)
# assert(delayTR>=0);

# delay_event_TR = pp.make_delay(delayTR)

# -----------------------
# Build sequence blocks
# -----------------------
# dummy loop
# for i in range(N_dummy_rep):

phase_areas = np.linspace(-Ny*delta_k/2, Ny*delta_k/2, Ny)
phase_areas_norm = phase_areas/np.max(phase_areas)

# adc = pp.make_adc(num_samples=Nx, 
#                   duration=gx_readout.flat_time, 
#                   delay=gx_readout.rise_time, 
#                   # phase_offset = 0,
#                   system=system)


adc_pos = pp.make_adc(num_samples=Nx/2, 
                  duration=gx_readout.flat_time, 
                  delay=gx_readout.rise_time, 
                  phase_offset = 0,
                  system=system)

adc_neg = pp.make_adc(num_samples=Nx/2, 
                  duration=gx_readout.flat_time, 
                  delay=gx_readout.rise_time, 
                  phase_offset = np.pi,
                  system=system)

# acquisition loop

### half alpha ###


gy_construct = pp.make_trapezoid(channel='y', 
                       area=np.max(phase_areas), 
                       #flat_time=3.2e-3, 
                       system=system)

# calculate TR:
blk_dur_01 = pp.calc_duration(rf_pos, gz_slicesel_pos)
blk_dur_02 = pp.calc_duration(gz_neg, gy_construct, gx_neg)
blk_dur_03 = pp.calc_duration(gx_readout, adc_pos)
TR = blk_dur_01 + 2*blk_dur_02 + blk_dur_03


seq.add_block(rf_excitation, gz_slicesel_pos)
blk_dur_init = pp.calc_duration(rf_excitation, gz_slicesel_pos)
delay = TR/2 - blk_dur_init
seq.add_block(pp.make_delay(delay))

for i in range(Ny):
  if i % 2 == 0:
    adc = adc_neg
    seq.add_block(rf_pos, gz_slicesel_pos)
    
  else:
    adc = adc_pos
    seq.add_block(rf_neg, gz_slicesel_pos)

  gy = pp.scale_grad(gy_construct, phase_areas_norm[i])
  # print(gy.rise_time + gy.flat_time + gy.fall_time)
  
  g1, g2, g3 =pp.align(left=[gz_neg, gy], right=[gx_neg])
  seq.add_block(g1,g2,g3)
  seq.add_block(gx_readout, adc)
  g1, g2, g3 =pp.align(left=[gx_neg, pp.scale_grad(gy, -1)], right=[gz_neg])
  #g1,g2,g3 =pp.align(left=[gx_neg, pp.scale_grad(pp.scale_grad(gy_construct, phase_areas_norm[i]), -1)], right=[gz_neg])
  seq.add_block(g1,g2,g3)
  #seq.add_block(gz_neg, gx_neg, pp.scale_grad(gy, -1))

  # seq.add_block(rf)
  # seq.add_block(adc, delay_event_TR)

# -----------------------
# Plot, timing check, write sequence
# -----------------------

# seq.plot(time_range=[0, 0.01])
# seq.plot()
seq.plot(show_blocks = True, time_range=[0, 0.02])
ok, error_report = seq.check_timing()
if ok:
    print("Timing check passed successfully.")
else:
    print("Timing check failed! Error listing follows:")
    for err in error_report:
        print(" -", err)

# Metadata
seq.set_definition(key='Name', value=name)

# Write to file
seq.write(name + '.seq')
print("Sequence written to " + name + ".seq")

# Optional: test report (not always implemented)
# try:
#     rep = seq.test_report()
#     for line in rep:
#         print(line), gy
# except AttributeError:
#     print("ℹ️ No test report available in this version of PyPulseq.")
