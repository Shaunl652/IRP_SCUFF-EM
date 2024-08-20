#!/usr/bin/env python3
import numpy as np
pi = np.pi
cos, sin, exp = np.cos, np.sin, np.exp
sqrt = np.sqrt
i = 1j

# Wavenumber and prefactor
wavelength = 1.55 
kmag = 2*pi/wavelength # Wave number
NumApp = 0.041 # numerical apperture, CHANGE THIS ONE
f0 = 10 # Filling factor from Maurer
w0 = wavelength/(pi*NumApp) # Beam waist
# These next two I've got no idea
# Think the first one is supposed to be the focal length so I've used that
Rinf = w0/(f0*NumApp) #1000.0*wavelength # is this prefactor even necessary for the numerics?
#prefactor = i*kmag*Rinf*exp(-i*kmag*Rinf)/(2*pi)


# Load points and weightings
# https://people.sc.fsu.edu/~jburkardt/datasets/sphere_lebedev_rule/sphere_lebedev_rule.html
txt = open("lebdev/lebedev_131.txt").read()
data = np.array([[float(v) for v in l.split()] for l in txt.strip().split("\n")])

# Sets up the bounds on theta and phi for the integration
def bounds(NA):
    phi = []
    theta = []
    w = []
    # From Maurer Appendix B, integration domain  uses theta_NA as the limits
    lim = np.rad2deg(np.arcsin(NA))
    for vals in data:
        if vals[1]<=lim: # vals[1] is the theta val?
            
            # Set up the phi values, should be 0->2pi
            phi.append(-(vals[0]-180)*pi/180) # shifts from -180...+180 to 0..360 while preserving sin & cos values
            # Set up thetas, does 0->theta_NA (shouldn't there be 2 bounds? (See App B))
            theta.append(vals[1]*pi/180)
            w.append(vals[2]) # Weight value of the angle
    return np.array(phi),np.array(theta),np.array(w)

phi,theta,w = bounds(NumApp)
# phi = -(data[:,0]-180)*pi/180 # shifts from -180...+180 to 0..360 while preserving sin & cos values
# theta = data[:,1]*pi/180
# w = data[:,2]

# normal vectors for each angle; see PNO Fig3.7 p59; Eq3.41 & Eq3.42

# 3.42
ntheta = np.array([cos(theta)*cos(phi),
                   cos(theta)*sin(phi),
                  -sin(theta)])

# 3.41
nphi   = np.array([ -sin(phi),
                    +cos(phi),
                   0*cos(phi)])

# Not sure
nrad = np.array([sin(theta)*cos(phi),
                  sin(theta)*sin(phi),
                  cos(theta)])

# # 3.40
# nrad = np.array([cos(phi),
#                  sin(phi),
#                  0*cos(theta)])

# https://en.wikipedia.org/wiki/Spherical_coordinate_system#Integration_and_differentiation_in_spherical_coordinates
Rmat = np.array([[sin(theta)*cos(phi), cos(theta)*cos(phi),-sin(phi)],
                 [sin(theta)*sin(phi), cos(theta)*sin(phi), cos(phi)],
                 [cos(theta)         ,-sin(theta),          0*phi  ]])

Rmat = np.array([                     [cos(theta)*cos(phi),-sin(phi)],
                                      [cos(theta)*sin(phi), cos(phi)],
                                      [-sin(theta),          0*phi  ]])

kvec = kmag*nrad

def make_beta(xvec):
    """Given a list of positions, xvec[3,N], return beta which can be used with Etilde to compute Efoc
        Efoc = np.einsum('abij,bj->ai', beta, Etilde)"""
    alpha = np.exp(1j*np.einsum('ai,aj->ij',xvec,kvec))
    beta = 4*pi*np.einsum('ij,j,abj->abij',alpha,w,Rmat)
    return beta
    

def Gaussian(w):
    """Make the Einf components, suitable for passing to E(Einf,rvecs), for a Gaussian beam of radius w"""
    # Arrays to store the field components
    N = len(theta)
    E = np.zeros((2,N),dtype=np.complex128)

    # Populate with a Gaussian beam
    Einc = lambda theta,phi : exp(-(Rinf*sin(theta)/w)**2)   # Eq3.52

    # Eq3.51
    E[0,:] = +Einc(theta,phi)*cos(phi)*sqrt(cos(theta)) # theta polarisation
    E[1,:] = -Einc(theta,phi)*sin(phi)*sqrt(cos(theta)) # phi polariation
    E[np.isnan(E)] = 0

    return E

def raster(Etilde):
    # Create a regular grid of points, as a 1D list
    import itertools as it
    x = np.linspace(-5, 5,100)*wavelength
    z = np.linspace(-5, 5,100)*wavelength
    xvec = np.array([[xi,0,zj] for xi in x for zj in z]).T

    # Find the beta value for this list of locations
    beta = make_beta(xvec)

    # Evaluate the field near to the focus
    Eval = np.einsum('abij,bj->ai', beta, Etilde)

    return x,z,Eval

def make_scuff(Einc):
    out = ["ASM"]
    for n in range(Einc.shape[1]):
        # Convert from polar to carteesian
        # E_theta * n_theta + E_phi * n_phi = [Ex,Ey,Ez]
        Ex,Ey,Ez = w[n]*(Einc[0,n]*ntheta[:,n] + Einc[1,n]*nphi[:,n])
        if abs(Ex)**2 + abs(Ey)**2 + abs(Ez)**2 == 0:
            continue
        # Builds the propogation vector
        nx,ny,nz = nrad[:,n]
        # Having these TWO lines means that we have counter propogating beams
        # Comment out the second line to get beam in positive z-direction
        out.append(f"   PW   {+nx:+.8f} {+ny:+.8f} {+nz:+.8f}   {np.real(Ex):+.8f}{np.imag(Ex):+.8f}i {np.real(Ey):+.8f}{np.imag(Ey):+.8f}i {np.real(Ez):+.8f}{np.imag(Ez):+.8f}i")
        out.append(f"   PW   {-nx:+.8f} {-ny:+.8f} {-nz:+.8f}   {np.real(Ex):+.8f}{np.imag(Ex):+.8f}i {np.real(Ey):+.8f}{np.imag(Ey):+.8f}i {np.real(Ez):+.8f}{np.imag(Ez):+.8f}i")
    out.append("END")
    return "\n".join(out)

if __name__=='__main__':

    # Compute Etilde for a Gaussian beam
    Etilde = Gaussian(w0)

    x,z,Eval = raster(Etilde)
    
    # Compute measures of this field
    Ival = np.einsum('ai,ai->i',Eval.conj(),Eval).real # Intensity
    Is = Ival.reshape((len(x),len(z)))                 # reshaped as the 2D grid

    Phase = np.angle(Eval[0,:]).reshape((len(x),len(z))) # Phase of the x polarisation component along x=0

    Gouy = np.unwrap(Phase[len(x)//2,:]) - kmag*z        # Phase with the kz dependence removed
    Gouy -= Gouy.mean()
    
    zR = np.pi*w0**2/1.55
    wz = w0*np.sqrt(1+(z/zR)**2)
    # Plot the field
    import matplotlib.pyplot as plt
    plt.figure()
    plt.title("Intensity")
    plt.pcolormesh(z,x,Is)
    plt.plot(x,0+wz,color='r',label='Paraxial')
    plt.plot(x,0-wz,color='r')
    plt.axis('equal')
    plt.legend()

    plt.figure()
    plt.title("Gouy phase")
    plt.plot(z,Gouy,label='Lebdev')
    plt.plot(z,np.arctan(z/zR),color='r',label='Paraxial')
    plt.legend()

    Ex = Eval[0,:].reshape((len(x),len(z)))
    plt.figure()
    plt.title(r"Re[$E_x$]")
    plt.pcolormesh(z,x,Ex.real)
    plt.plot(x,0+wz,color='r',label='Paraxial')
    plt.plot(x,0-wz,color='r')
    plt.legend()
    
    plt.axis('equal')    
    
    plt.show()

    output = make_scuff(Etilde)
    
    with open('sources_ASM.list','w') as f:
        f.write(output)
