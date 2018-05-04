function gf = relaxEpithelium2D(gi, maxIterations)
    % relax input lattice to mechanical equilibrium
    %
    % g = relaxEpithelium(ginit)
    %
    % wrapper for mex code to do some easy input checking
    % this stuff is messy in C and makes the code less readable
    
    if nargin == 1
        maxIterations = 10;
    end
    
    assert(isstruct(gi), 'input needs to be a structure');
    
    % check lattice structure has the right fields
    assert(isfield(gi,'cells'), 'ginit needs field cells');
    assert(isfield(gi,'bonds'), 'ginit needs field bonds');
    assert(isfield(gi,'verts'), 'ginit needs field verts');
    
    % check fields are of right type
    assert(iscell(gi.cells), 'ginit.cells should be cell array');
    assert(isnumeric(gi.bonds) && size(gi.bonds,2)==4,...
                                    'ginit.bonds should be Nx4 matrix');
    assert(isnumeric(gi.verts) && size(gi.verts,2)==3,...
                                    'ginit.verts should be Nx3 matrix');	
    
    [v, status, E, tForce, pForce, perimForce, A, T, cellForce] =...
                            Lattmin(    gi.verts, gi.bonds, gi.cells,...
                                        gi.param,...
                                        gi.T0, gi.A0, gi.p0, gi.l0,...
                                        gi.pT0, gi.kA0, maxIterations);

	% update lattice
    gf = gi;
    gf.verts = v;
    
    % output useful for debugging
    gf.tForce = tForce;
    gf.pForce = pForce;
    gf.perimForce = perimForce;
    gf.A = A;
    gf.T = T;
    gf.E = E;
    gf.cellForce = cellForce;
    
    %gf.stress = computeDilationStress2(gf);
end