import meep as mp
import numpy as np
import argparse

def main(args):
    resolution = 50         # pixels/um

    normalize = args.normalize
    particle = args.particle

    r = args.r              # hole radius
    w_start = args.w_start  # center width
    w_end = args.w_end      # wvg width
    h = args.hh             # height of wvg                               
    a = args.a              # period of holes
    r_par = args.r_par      # radius of particle
    offset = args.offset    # offset of the particle center from cavity center

    pad = 1                 # additional padding at ends of evg
    dair = 1                # air region
    dpml = 1                # PML thickness
    pml_layers = [mp.PML(dpml)]
    Ndef = args.Ndef        # number of defects
    Nwvg = args.Nwvg        # number of wvg cells

    sx = dpml + pad + 2*(Ndef+Nwvg)*a + pad + dpml
    sy = dpml+dair+w_start+dair+dpml
    sz = dpml+dair+h+dair+dpml
    
    cell = mp.Vector3(sx,sy,sz)

    def taper(i):
        return w_start + i*i*(w_end-w_start)/(Ndef*Ndef)  # cubic taper of cell width

    # materials
    Si = mp.Medium(index=3.46)      # @1550 nm
    SiO2 = mp.Medium(index=1.44)    # @1550 nm

    # construct the geometry objects
    geometry = []
    # input wvg    
    geometry.append(mp.Block(material=Si, center=mp.Vector3(), size=mp.Vector3(mp.inf,w_end,h)))
    
    # if it is not normalization run, make nanobeam cavity
    if not normalize: 
        for mm in range(Ndef):
            geometry.append(mp.Block(material=Si, center=mp.Vector3(), size=mp.Vector3((2*mm+1)*a,taper(mm),h)))    # taper section
    
        for mm in range(Ndef+Nwvg):
            geometry.append(mp.Cylinder(material=mp.air, radius=r, height=mp.inf, center=mp.Vector3(+mm*a,0,0)))
            geometry.append(mp.Cylinder(material=mp.air, radius=r, height=mp.inf, center=mp.Vector3(-mm*a,0,0)))    # taper/mirror section
        
    if particle:
            p_center = mp.Vector3(0.00, offset, 0.00) # place the particle, where you want
            geometry.append(mp.Sphere(material=SiO2, radius=r_par, center=p_center))                                # add particle

    # source
    lam_min = 1.45
    lam_max = 1.65
    fmin = 1/lam_max
    fmax = 1/lam_min
    fcen = 0.5*(fmin + fmax)
    df = fmax-fmin

    symmetries = []
    
    sources = [mp.Source(mp.GaussianSource(fcen, fwidth=df),
                                     component=mp.Ey,
                                     center=mp.Vector3(-0.5*sx+dpml),
                                     size=mp.Vector3(0,w_end,h))]

    # construct sim class
    sim = mp.Simulation(resolution=resolution,
                        cell_size=cell,
                        boundary_layers=pml_layers,
                        sources=sources,
                        geometry=geometry,
                        symmetries=symmetries)

    # accumulate fft fields at the end of output wvg    
    freg = mp.FluxRegion(center=mp.Vector3(0.5*sx-dpml-pad/2), size=mp.Vector3(0,2*w_end,2*h))
    lam_min = 1.545
    lam_max = 1.560
    fmin = 1/lam_max
    fmax = 1/lam_min
    fcen = 0.5*(fmin + fmax)
    df = fmax-fmin
    nfreq = 500 # number of frequencies at which to compute flux

    # add transmitted flux
    trans = sim.add_flux(fcen, df, nfreq, freg)

    if particle:
            r_offset = r_par + 2/resolution # offset btw particle and forcebox ~ 1px
            N_planes = 9
            px_plane = 2*r_offset/N_planes 
                        
            force_reg1  = mp.ForceRegion(center=p_center + mp.Vector3(x=+r_offset                          ), direction=mp.X, weight=+1, size=mp.Vector3(y=5*px_plane, z=2*r_offset))
            force_reg2  = mp.ForceRegion(center=p_center + mp.Vector3(x=+r_offset-px_plane,   y=-3*px_plane), direction=mp.X, weight=+1, size=mp.Vector3(y=px_plane,   z=2*r_offset))
            force_reg3  = mp.ForceRegion(center=p_center + mp.Vector3(x=+r_offset-2*px_plane, y=-4*px_plane), direction=mp.X, weight=+1, size=mp.Vector3(y=px_plane,   z=2*r_offset))
            force_reg4  = mp.ForceRegion(center=p_center + mp.Vector3(x=+r_offset-px_plane,   y=+3*px_plane), direction=mp.X, weight=+1, size=mp.Vector3(y=px_plane,   z=2*r_offset))
            force_reg5  = mp.ForceRegion(center=p_center + mp.Vector3(x=+r_offset-2*px_plane, y=+4*px_plane), direction=mp.X, weight=+1, size=mp.Vector3(y=px_plane,   z=2*r_offset))
            force_reg6  = mp.ForceRegion(center=p_center + mp.Vector3(x=-r_offset+px_plane,   y=-3*px_plane), direction=mp.X, weight=-1, size=mp.Vector3(y=px_plane,   z=2*r_offset))
            force_reg7  = mp.ForceRegion(center=p_center + mp.Vector3(x=-r_offset+2*px_plane, y=-4*px_plane), direction=mp.X, weight=-1, size=mp.Vector3(y=px_plane,   z=2*r_offset))
            force_reg8  = mp.ForceRegion(center=p_center + mp.Vector3(x=-r_offset                          ), direction=mp.X, weight=-1, size=mp.Vector3(y=5*px_plane, z=2*r_offset))
            force_reg9  = mp.ForceRegion(center=p_center + mp.Vector3(x=-r_offset+px_plane,   y=+3*px_plane), direction=mp.X, weight=-1, size=mp.Vector3(y=px_plane,   z=2*r_offset))
            force_reg10 = mp.ForceRegion(center=p_center + mp.Vector3(x=-r_offset+2*px_plane, y=+4*px_plane), direction=mp.X, weight=-1, size=mp.Vector3(y=px_plane,   z=2*r_offset))
            
            force_reg11 = mp.ForceRegion(center=p_center + mp.Vector3(                          y=+r_offset), direction=mp.Y, weight=+1, size=mp.Vector3(x=5*px_plane, z=2*r_offset))
            force_reg12 = mp.ForceRegion(center=p_center + mp.Vector3(x=-3*px_plane,   y=+r_offset-px_plane), direction=mp.Y, weight=+1, size=mp.Vector3(x=px_plane,   z=2*r_offset))
            force_reg13 = mp.ForceRegion(center=p_center + mp.Vector3(x=-4*px_plane, y=+r_offset-2*px_plane), direction=mp.Y, weight=+1, size=mp.Vector3(x=px_plane,   z=2*r_offset))
            force_reg14 = mp.ForceRegion(center=p_center + mp.Vector3(x=+3*px_plane,   y=+r_offset-px_plane), direction=mp.Y, weight=+1, size=mp.Vector3(x=px_plane,   z=2*r_offset))
            force_reg15 = mp.ForceRegion(center=p_center + mp.Vector3(x=+4*px_plane, y=+r_offset-2*px_plane), direction=mp.Y, weight=+1, size=mp.Vector3(x=px_plane,   z=2*r_offset))
            force_reg16 = mp.ForceRegion(center=p_center + mp.Vector3(x=-3*px_plane,   y=-r_offset+px_plane), direction=mp.Y, weight=-1, size=mp.Vector3(x=px_plane,   z=2*r_offset))
            force_reg17 = mp.ForceRegion(center=p_center + mp.Vector3(x=-4*px_plane, y=-r_offset+2*px_plane), direction=mp.Y, weight=-1, size=mp.Vector3(x=px_plane,   z=2*r_offset))
            force_reg18 = mp.ForceRegion(center=p_center + mp.Vector3(                          y=-r_offset), direction=mp.Y, weight=-1, size=mp.Vector3(x=5*px_plane, z=2*r_offset))
            force_reg19 = mp.ForceRegion(center=p_center + mp.Vector3(x=+3*px_plane,   y=-r_offset+px_plane), direction=mp.Y, weight=-1, size=mp.Vector3(x=px_plane,   z=2*r_offset))
            force_reg20 = mp.ForceRegion(center=p_center + mp.Vector3(x=+4*px_plane, y=-r_offset+2*px_plane), direction=mp.Y, weight=-1, size=mp.Vector3(x=px_plane,   z=2*r_offset))

            force_reg21 = mp.ForceRegion(center=p_center + mp.Vector3(               z=+r_offset), direction=mp.Z, weight=+1, size=mp.Vector3(x = 5*px_plane, y=9*px_plane))
            force_reg22 = mp.ForceRegion(center=p_center + mp.Vector3(x=+3*px_plane, z=+r_offset), direction=mp.Z, weight=+1, size=mp.Vector3(x = px_plane,   y=7*px_plane))
            force_reg23 = mp.ForceRegion(center=p_center + mp.Vector3(x=-3*px_plane, z=+r_offset), direction=mp.Z, weight=+1, size=mp.Vector3(x = px_plane,   y=7*px_plane))
            force_reg24 = mp.ForceRegion(center=p_center + mp.Vector3(x=+4*px_plane, z=+r_offset), direction=mp.Z, weight=+1, size=mp.Vector3(x = px_plane,   y=5*px_plane))
            force_reg25 = mp.ForceRegion(center=p_center + mp.Vector3(x=-4*px_plane, z=+r_offset), direction=mp.Z, weight=+1, size=mp.Vector3(x = px_plane,   y=5*px_plane))
            force_reg26 = mp.ForceRegion(center=p_center + mp.Vector3(               z=-r_offset), direction=mp.Z, weight=-1, size=mp.Vector3(x = 5*px_plane, y=9*px_plane))
            force_reg27 = mp.ForceRegion(center=p_center + mp.Vector3(x=+3*px_plane, z=-r_offset), direction=mp.Z, weight=-1, size=mp.Vector3(x = px_plane,   y=7*px_plane))
            force_reg28 = mp.ForceRegion(center=p_center + mp.Vector3(x=-3*px_plane, z=-r_offset), direction=mp.Z, weight=-1, size=mp.Vector3(x = px_plane,   y=7*px_plane))
            force_reg29 = mp.ForceRegion(center=p_center + mp.Vector3(x=+4*px_plane, z=-r_offset), direction=mp.Z, weight=-1, size=mp.Vector3(x = px_plane,   y=5*px_plane))
            force_reg30 = mp.ForceRegion(center=p_center + mp.Vector3(x=-4*px_plane, z=-r_offset), direction=mp.Z, weight=-1, size=mp.Vector3(x = px_plane,   y=5*px_plane))
            
            par_force_x = sim.add_force(fcen, df, nfreq, force_reg1, force_reg2, force_reg3,force_reg4, force_reg5, force_reg6, force_reg7, force_reg8, force_reg9, force_reg10)
            par_force_y = sim.add_force(fcen, df, nfreq, force_reg11, force_reg12, force_reg13, force_reg14, force_reg15, force_reg16, force_reg17, force_reg18, force_reg19, force_reg20)
            par_force_z = sim.add_force(fcen, df, nfreq, force_reg21, force_reg22, force_reg23, force_reg24, force_reg25, force_reg26, force_reg27, force_reg28, force_reg29, force_reg30)  

    sim.use_output_directory()

    # output planes
    volXY = mp.Volume(center=mp.Vector3(), size=mp.Vector3(sx,sy,0))
    volXZ = mp.Volume(center=mp.Vector3(), size=mp.Vector3(sx,0,sz))

    sim.run(mp.in_volume(volXY, mp.with_prefix('XY', mp.at_beginning(mp.output_epsilon))),
            mp.in_volume(volXZ, mp.with_prefix('XZ', mp.at_beginning(mp.output_epsilon))),
            mp.in_volume(volXY, mp.with_prefix('XY', mp.at_every(1000, mp.synchronized_magnetic(mp.output_efield)))),
            mp.in_volume(volXZ, mp.with_prefix('XZ', mp.at_every(1000, mp.synchronized_magnetic(mp.output_efield)))),
            until_after_sources=mp.stop_when_fields_decayed(50, mp.Ey, mp.Vector3(0.5*sx-dpml-pad/2), 1e-6))
    
    # get the fluxes and frequencies
    tran_flux = mp.get_fluxes(trans)
    
    if particle:
        flux_freqs = mp.get_force_freqs(par_force_x)
        force_y = mp.get_forces(par_force_y)
        force_x = mp.get_forces(par_force_x)
        force_z = mp.get_forces(par_force_z)    
            
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-normalize', action='store_true', default=False, help='run the normalization? (default: False)')
    parser.add_argument('-particle', action='store_true', default=False, help='run with the particle? (default: False)')
    parser.add_argument('-offset', type=float, default=0.0, help='Offset of  the particle from cavity center')
    parser.add_argument('-w_start', type=float, default=0.7, help='starting width (default: 0.70 um)')
    parser.add_argument('-w_end', type=float, default=0.5, help='ending width (default: 0.50 um)')
    parser.add_argument('-a', type=float, default=0.52, help='period of PC (default: 0.52 um)')
    parser.add_argument('-r', type=float, default=0.15, help='hole radius (default: 0.15 um)')
    parser.add_argument('-r_par', type=float, default=0.05, help='radius of particle(default: 0.05 um)')
    parser.add_argument('-hh', type=float, default=0.22, help='waveguide height (default: 0.22 um)')
    parser.add_argument('-Ndef', type=int, default=7, help='number of defect periods (default: 7)')
    parser.add_argument('-Nwvg', type=int, default=0, help='number of defect periods (default: 0)')
    args = parser.parse_args()
    main(args)
    