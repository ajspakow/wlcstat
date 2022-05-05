# -*- coding: utf-8 -*-
r"""Wigner D class

This class allows you to quickly calculate Wigner D matrix values using the symbolic lookup Wigner D
functionality available from sympy.physics. It will be slow when first calculating a set of l, m, 
and j values but will hence-forth become fast and evaluate quickly for different angles.

This class also normalizes the standard Wigner D matrix values according to Andy's notes:

.. math::

    \begin{equation*}
    {\mathcal{D}}^{mj}_{l} = \sqrt{\frac{2l+1}{8\pi^2}}D^{mj}_{l}
    \end{equation*}

Here, :math:`{\mathcal{D}}^{mj}_{l}` refers to Andy's renormalized Wigner D's, and :math:`D^{mj}_{l}`
refers to the standard wigner D functions found in the sympy docs and on Wikipedia.
This normalization ensures the wigner D functions can be treated as normalized probability distributions 
that obey the following ladder operation condition:

.. math::

    \begin{equation*}
    \cos{\theta}{\mathcal{D}}^{mj}_{l} = \alpha^{mj}_{l}{\mathcal{D}}^{mj}_{l-1} +
    \beta^{mj}_{l}{\mathcal{D}}^{mj}_{l} + \alpha^{mj}_{l+1}{\mathcal{D}}^{mj}_{l+1}
    \end{equation*}

where :math:`\alpha^{mj}_{l}=\sqrt{\frac{(l-m)(l+m)(l-j)(l+j)}{l^2(4l^2-1)}}` and :math:`\beta^{mj}_{l}=\frac{mj}{l(l+1)}`
come from the Clebsh-Gordon coefficients.

Note that sympy docs uses :math:`j, m, m'` instead of :math:`l, m, j`. We will use the latter notation
to be consistent with Andy's paper on Wormlike chain statistics with twist and fixed ends.

Example
-------
.. code-block:: python

    >>> import wignerD as wd
    >>> myset = wd.wigner_d_vals()
    >>> ans = myset.get(l,m,j,alpha,beta,gamma)
    >>> mat = myset.get_mtrx(l,m,j,alpha,beta,gamma)

You can also use a faster algorithm if you have :math:`\gamma=0`.

.. code-block:: python

    >>> ans = myset.axial(l,m,j,alpha,beta)

This class builds a dictionary of functions, with keys (l, m, j) mapping to functions of (alpha, beta, gamma).
To use these functions directly (without the get() or axial() methods),

.. code-block:: python

    >>> l, m, j = 1, 0, 0
    >>> alpha, beta, gamma = 0, np.pi, 0
    >>> myset.fun[(l, m, j)](alpha, beta, gamma)

"""

from sympy.physics.quantum.spin import Rotation
from sympy import symbols
import numpy as np
from sympy.utilities import lambdify

class wigner_d_vals:
    def __init__(self):
        self.fun = {}
        """dictionary of functions of (alpha, beta, gamma), key: (l,m,j)"""
        self.fun_axial = {}
        """dictionary of functions, key: (l,m,j) where :math:`\gamma=0`; axial symmetry"""
    
    def normalize(self,l,m,j):
        """Returns the normalization constant :math:`\sqrt{\frac{2l+1}{8\pi^2}}`"""
        return np.sqrt((2*l+1)/(8*np.pi**2))

    def axial(self,l,m,j,alpha=None,beta=None,norm=True):
        """Returns a Wigner D matrix element for a defined rotation given by the Euler angles
        :math:`\alpha, \beta, \gamma = 0.` If alpha, beta are not specified, simply updates the dictionary
        self.fun_axial to include a function of (alpha, beta) for the key (l, m, j).

        Parameters
        ----------
        l : int
            total angular momentum
        m : int
            eigenvalue of angular momentum along axis after rotation
        j : int
            eigenvalue of angular momentum along rotated axis
        alpha : float, symbol
            the final rotation about z''
        beta : float, symbol
            the second rotation about y'
        gamma : float, symbol
            the initial rotation about z
        norm: bool (default = True)
            whether or not to normalize to match Andy's wigner D's

        """
        key = (l,m,j)
        if abs(m)>l or abs(j)>l:
            raise Exception('m must be <= l')
        if key in self.fun_axial:
            fun=self.fun_axial[key]
        else: #update self.fun_axial dictionary to add key (l, m, j)
            al, be = symbols('alpha beta')
            fun=lambdify((al,be),Rotation.D(l,m,j,al,be,0).doit(),"numpy")
            self.fun_axial[key]=fun
        if (alpha is None) != (beta is None):
            raise ValueError('Either none or both of alpha and beta must be specified.')
        #only return a value if Euler angles were specified
        if alpha is not None:
            if norm:
                return self.normalize(l,m,j)*fun(alpha, beta)
            else:
                return fun(alpha,beta)    
            
    def get(self,l,m,j,alpha=None,beta=None,gamma=None,norm=True):
        """Returns a Wigner D matrix element for a defined rotation given by the Euler angles
        :math:`\alpha, \beta, \gamma.` If alpha, beta, gamma are not specified, simply updates the dictionary
        self.fun to include a function of (alpha, beta, gamma) for the key (l, m, j). 

        Parameters
        ----------
        l : int
            total angular momentum
        m : int
            eigenvalue of angular momentum along axis after rotation
        j : int
            eigenvalue of angular momentum along rotated axis
        alpha : float, symbol
            the final rotation about z''
        beta : float, symbol
            the second rotation about y'
        gamma : float, symbol
            the initial rotation about z
        norm: bool (default = True)
            whether or not to normalize to match Andy's wigner D's

        """
        key = (l,m,j)
        if abs(m)>l or abs(j)>l:
            raise Exception('m must be <= l')
        if key in self.fun:
            fun=self.fun[key]
        else: #update self.fun dictionary to add key (l, m, j)
            al, be, ga = symbols('alpha beta gamma')
            fun=lambdify((al,be,ga),Rotation.D(l,m,j,al,be,ga).doit(),"numpy")
            self.fun[key]=fun
        if (alpha is None) != (beta is None) != (gamma is None):
            raise ValueError('Either none or all of alpha, beta, gamma must be specified.')
        #only return a value if Euler angles were specified
        if alpha is not None:
            if norm:
                return self.normalize(l,m,j)*fun(alpha, beta, gamma)
            else: 
                return fun(alpha,beta,gamma)             

    def get_mtrx(self,alpha,beta,gamma,nlam,norm=True):
        """Returns a list representation of a Wigner D matrix for a defined rotation given by the Euler angles
        :math:`\alpha, \beta, \gamma.`

        Parameters
        ----------
        l : int
            total angular momentum
        m : int
            eigenvalue of angular momentum along axis after rotation
        j : int
            eigenvalue of angular momentum along rotated axis
        alpha : float
            the final rotation about z''
        beta : float
            the second rotation about y'
        gamma : float
            the initial rotation about z
        nlam : int
            number of l values to compute (maximum l = nlam-1)
        norm: bool (default = True)
            whether or not to normalize to match Andy's wigner D's
        
        Returns
        -------
        Dlist : (nlam,) list
            List indexed by l, such that Dlist[l] is a (2l+1) by (2l+1) matrix
            indexed by [m-l, j-l] where m,j range from -l to l
        """
        out=[] # list of matricies [l][m-l,j-l]
        for l in range(0,nlam):
            print('l',l)
            D = np.zeros((2*l+1,2*l+1),dtype='complex')
            for j in range(-l,l+1):
                for m in range(-l,l+1):
                    D[m-l,j-l]=self.get(l,m,j,alpha,beta,gamma, norm)
            out.append(D)
        return out

