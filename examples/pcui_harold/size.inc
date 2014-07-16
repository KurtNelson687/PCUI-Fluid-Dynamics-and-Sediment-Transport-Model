C...... ni/px must be dividable by 16
C...... nj/py must be dividable by 16
C...... nk/pz must be dividable by 16

        implicit none

	integer, parameter :: px = 64
	integer, parameter :: py = 4
	integer, parameter :: pz = 3 

	integer, parameter :: ni = 2048
	integer, parameter :: nj = 128 
	integer, parameter :: nk = 96

	integer, parameter :: ni1 = ni / 2
	integer, parameter :: nj1 = nj / 2
	integer, parameter :: nk1 = nk / 2
	integer, parameter :: ni2 = ni / 4
	integer, parameter :: nj2 = nj / 4
	integer, parameter :: nk2 = nk / 4
	integer, parameter :: ni3 = ni / 8
	integer, parameter :: nj3 = nj / 8
	integer, parameter :: nk3 = nk / 8
	integer, parameter :: ni4 = ni / 16
	integer, parameter :: nj4 = nj / 16
	integer, parameter :: nk4 = nk / 16

	integer, parameter :: nni = ni / px
	integer, parameter :: nnj = nj / py
	integer, parameter :: nnk = nk / pz
	integer, parameter :: nni1 = ni1 / px
	integer, parameter :: nnj1 = nj1 / py
	integer, parameter :: nnk1 = nk1 / pz
	integer, parameter :: nni2 = ni2 / px
	integer, parameter :: nnj2 = nj2 / py
	integer, parameter :: nnk2 = nk2 / pz
	integer, parameter :: nni3 = ni3 / px
	integer, parameter :: nnj3 = nj3 / py
	integer, parameter :: nnk3 = nk3 / pz
	integer, parameter :: nni4 = ni4 / px
	integer, parameter :: nnj4 = nj4 / py
	integer, parameter :: nnk4 = nk4 / pz
