	double precision bx, by, bz, dm, cm, am, bm,
     <    maxRatio, numWallUnits
	integer stretchx, stretchy, stretchz
	integer stopStretch
	double precision wallUnit, stretchFactor,
     <		dx, dz, dy, dyMin, domainVol
	double precision, dimension(nj) ::
     <		dyAll, yAll
	common/mapping/
     <		bx, by, bz, dm, cm, am, bm,dyAll,
     <		stretchx, stretchy, stretchz, yAll,
     <          domainVol

