program testTimestreamInfo
use healpix_types
use pix_tools
use fitstools
use quiet_utils
use l2_fileutils
use quiet_fileutils
use quiet_calib_mod
use quiet_noise_mod
use ds_types
use ds_azimuth
use ds_oslo_fitstools
use ds_multidetector
use ds_solver
use ds_fitstools

type(l2_filelist),   pointer, dimension(:) :: target_list
integer(i4b), pointer, dimension(:) :: module_list, module_list_cut
integer unit
integer t
integer nmodules
type(ds_timestreamInfo) TS

unit=41
call initialize_quiet_db_mod(unit, "runlist.txt", "/data4/quiet/joez/test1", "/data4/quiet/joez/")
call read_module_list(unit, "auxilliary/module_list.txt", module_list)
call initialize_accepted_scan_list(unit,"accepted_4a.txt")
call get_target_list(1,"cmb", "patch_gc", target_list)
call removeTemperatureFromModuleList(module_list, module_list_cut, nmodules)
call buildTimestreamInfo("accepted_4a.txt",target_list,TS,module_list_cut)
do t=0,TS%n-1
	write(*,*) TS%run(t), TS%scan(t), TS%mod(t), TS%pair(t), TS%needs_azfilter(t), trim(TS%filename(t))
enddo
end program testTimestreamInfo