	double precision bx, by, bz, dm, cm, am, bm
	integer stretchx, stretchy, stretchz
	integer maxRatio, numWallUnits, stopStretch
	double precision wallUnit, stretchFactor,
     <		dx, dz, dy, dyMin
	double precision, dimension(nj) ::
     <		dyAll, yAll
	common/mapping/
     <		bx, by, bz, dm, cm, am, bm,
     <		stretchx, stretchy, stretchz, yAll

