# This is a python program to generate the options setting and examples
# in the destriper.
#
# It will certainly work as far back as python2.5, and possibly 2.4.
#

import StringIO

DEFAULT_PARAMETER_FILENAME = "default_params.ini"
FORTRAN_FILENAME = "ds_cbass_options.f90"

#(Parameter_name,Parameter_type,Default_value,Help_text)

basic_options=[
("file_list",str,"files.txt","List of FITS files to make maps from"),
("nside",int,128,"The Healpix NSIDE parameter of the output maps"),
("offsetLength",int,1000,"The length of the destriping offsets"),
("naive_mode",bool,False,"Make a naive map without destriping"),
("traditional_mode",bool,True,"Do not use noise information when destriping"),
("do_temperature",bool,True,"Make temperature maps"),
("do_polarization",bool,False,"Make temperature maps"),
("make_azel_maps",bool,False,"Make maps in azimuth/elevation coordinates instead of ra/dec"),
("cut_extension_name",str,"","Specify the name of the FITS extension in the files that contains the cut information.  Default is to use the first one"),
("polarization_angle",float,0.0,"The polarization offset angle of the telescope."),
("enable_I1",bool,True,"The I1 channel is included in the data (if do_temperature = T)"),
("enable_I2",bool,True,"The I2 channel is included in the data (if do_temperature = T)"),
("enable_Q1",bool,True,"The Q1 channel is included in the data (if do_polarization = T)"),
("enable_Q2",bool,True,"The Q2 channel is included in the data (if do_polarization = T)"),
("enable_U1",bool,True,"The U1 channel is included in the data (if do_polarization = T)"),
("enable_U2",bool,True,"The U2 channel is included in the data (if do_polarization = T)"),
]

output_options=[
("verbosity",int,2,"Verbosity level of code.  0=silent, 1=quiet, 2=noisy, 3=debug.  Default is 2."),
("output_filename",str,"map.fits","Filename for the output maps"),
("save_hits",bool,True,"Save a hit count map"),
("hits_filename",str,"map.hits","Filename for the output hit count map"),
("save_covariance",bool,True,"Save the covariance of the output map"),
("covariance_filename",str,"map.cov","Filename for the output covariance"),
("save_offsets",bool,False,"Save the offsets fitted to the timestream in destriping"),
("offset_dir",str,"offsets/","Directory in which to save the offsets"),
("ring_ordering",bool,True,"Whether to use Healpix ring ordering for the output map"),
("save_name_mapping",bool,False,"Whether to save the mapping of file names to ID numbers in a text file"),
("name_mapping_file",str,"name_mappings.txt","Filename for the mapping of file names to ID numbers"),
]

tangent_plane_mapping=[
('plane_map',bool,False,"Make a map of a small sky section by projecting onto a tangent plane."),
('plane_center_ra',float,0.0,'RA of center of plane map'),
('plane_center_dec',float,0.0,'DEC of center of plane map'),
('plane_pixel_size',float,0.1,'Pixel size in degrees for plane map'),
('plane_size',int,100,'Length of one side of plane map in pixels'),
]

timestream_options=[
("should_cut_timestreams",bool,False,"Cut down the lengths of the timestreams, for debugging."),
("cutfraction",float,1.0,"The fraction of the timestreams to remove, if should_cut_timestreams is true"),
("simulate",bool,False,"Replace the loaded timestreams with simulations"),
("simulation_filename",str,"simulated_map.fits","The signal map to get the simulated signal from"),
("subtractMeans",bool,True,"Subtract the mean from each timestream chunk before processing"),
("subtract_linear",bool,False,"Subtract a linear trend from each timestream before processing"),
("use_tod_filters",bool,False,"Filter the TOD before processing (Not functional now)"),
("partial_out",bool,False,"Output partial healpix maps with just the hit pixels tabulated (cuts down file size for maps of small sky regions)"),
]

pcg_options=[
("pcg_tolerance",float,1.0e-6,"The tolerance of the PCG solver."),
("data_prior",bool,False,"Construct the prior on the timestream noise directly from the data"),
("do_sigma_cut",bool,False,"Remove troublesome scans which are preventing convergence of the destriper, usuall because they give very different results to other scans of the same pixels"),
("sigma_cut_value",float,100.0,"The number of sigmas away from the rest a scan must be to be considered 'troublesome'"),
("sigma_cut_time",int,10,"The number of iterations at which to cut remove troublesome scans"),
]

azimuth_options=[
("use_azimuth_offsets",bool,False,"Attempt to remove azimuth-correlated noise (Not functioning)"),
("azimuth_nbins",int,20,"The number of bins in azimuth to use for azimuth correlated noise (Not functional)"),
("azimuth_corr_time",float,20.0,"The maximum correlation time of azimuth noise"),
("simulate_az",bool,False,"Add simulated azimuth-correlated noise (Not functional)"),
]

options=[
("Basic Options",basic_options),
("Ouput Control Options",output_options),
("Plane Mapping",tangent_plane_mapping),
("Timestream Processing Options",timestream_options),
("Destriping Conjugate Gradient Options",pcg_options),
("Azimuth Noise Options (Not functional)",azimuth_options),
]



def fortran_string(value,ptype):
	if ptype==str:
		return '"%s"'%value
	elif ptype in (float,int):
		return str(value)
	elif ptype==bool:
		if value:
			return '.true.'
		else:
			return '.false.'
	else:
		raise ValueError("Unknown parameter type")

def declaration_string(ptype):
	if ptype==str:
		return 'character(256) :: '
	elif ptype == float:
		return 'real(dp) :: '
	elif ptype == int:
		return 'integer :: '
	elif ptype==bool:
		return 'logical :: '
	else:
		raise ValueError("Unknown parameter type")
	
FITS_CODES = {
		float : "D",
		int : "J",
		bool  : "L",
		str : "S"
	}

FORTRAN_CODES={
	float:'E10.5',
	str:'A',
	int:'I0',
	bool:'L'
}	

def ini_string(value,ptype):
	if ptype in (float,int,str):
		return str(value)
	elif ptype==bool:
		if value:
			return 'T'
		else:
			return 'F'
	else:
		raise ValueError("Unknown parameter type")


def ini_call_for_type(ptype):
	if ptype==str:
		return "ini_read_string_file2"
	elif ptype==float:
		return "ini_read_double_file"
	elif ptype==int:
		return "ini_read_int_file"
	elif ptype==bool:
		return "ini_read_logical_file"
	else:
		raise ValueError("Unknown parameter type")
	

def generate_fortran(options):
	S=StringIO.StringIO()
	S.write("""! DO NOT EDIT THIS FILE - IT IS GENERATED FROM build_options.py. Edit that.

module ds_cbass_option_utils
use inifile
use ds_precision
use ds_utils
implicit none

type ds_cbass_options
""")

	for set_name,option_set in options:
		S.write("\n\t!!!!!!!!!!!!!!!!!\n")		
		S.write("\t!!! %s\n"%set_name)		
		S.write("\t!!!!!!!!!!!!!!!!!\n\n")
		for name,ptype,default,help in option_set:
			S.write("\t!%s\n"%help)
			S.write("\t%s%s\n\n"%(declaration_string(ptype),name))
	S.write("end type ds_cbass_options\n\n\n")
	S.write("contains\n\n")

	S.write("subroutine set_default_options(options)\n")
	S.write("\ttype(ds_cbass_options) :: options\n")
	for set_name,option_set in options:
		for name,ptype,default,help in option_set:
			S.write("\toptions%%%s = %s\n"%(name,fortran_string(default,ptype)) )
	S.write("end subroutine set_default_options\n\n\n")
	
	S.write("""
subroutine read_options(filename,options)
	type(ds_cbass_options) :: options
	character(*) :: filename
	Type(TIniFile) :: Ini
	integer unit
	logical error
	
	unit=ds_get_lun()
	call Ini_Open_File(ini, trim(filename), unit,  error)
  	if (error) then
		write(*,*) 'Unable to open or read file: '//trim(filename)
		stop
	endif
	
	call set_default_options(options)

""")
	for set_name,option_set in options:
		for name,ptype,default,help in option_set:
			S.write("\toptions%%%s = %s(ini,'%s',options%%%s)\n" % (name,ini_call_for_type(ptype),name,name ) )
	S.write("end subroutine read_options\n\n\n")
	
	S.write("""
subroutine options_to_fits_cards(options,cards)
	use fl_lists
	type(ds_cbass_options) :: options
	integer unit
	type(fl_string80_list) :: cards
	character(256) :: template
	character(80) :: card
	integer status, hdutype
	status=0
	hdutype=0
	template='COMMENT  ------------------------'
	call FTGTHD(trim(template),card,hdutype,status)
	if (status .ne. 0) then 
		write(*,'(A,I0,A)') "ERROR ",status," PARSING:"
		write(*,'(A)') trim(template)
	endif
	call fl_append_string80_list(cards,card)
	template='COMMENT  ---LIST OF MAPMAKING PARAMETERS---'
	call FTGTHD(trim(template),card,hdutype,status)
	if (status .ne. 0) then 
		write(*,'(A,I0,A)') "ERROR ",status," PARSING:"
		write(*,'(A)') trim(template)
	endif
	call fl_append_string80_list(cards,card)
	template='COMMENT  ------------------------'
	call FTGTHD(trim(template),card,hdutype,status)
	if (status .ne. 0) then 
		write(*,'(A,I0,A)') "ERROR ",status," PARSING:"
		write(*,'(A)') trim(template)
	endif
	call fl_append_string80_list(cards,card)
	
""")
	i=0
	for set_name,option_set in options:
		for (name,ptype,default,help) in option_set:
			fits_card(S,i,name,ptype,help)
			i+=1
	S.write("end subroutine options_to_fits_cards\n")		
	S.write("end module ds_cbass_option_utils\n")
	S.seek(0)
	return S.read()





def fits_card(S,i,name,ptype,help):
	key_keyword='DSKEY%.3d'%i
	key_comment=fortran_string(' / Descart parameter %d name' % i,str)
	key_value=fortran_string(name,str)
	
	template_write = "	write(template,'(A,A,A)') '%s  ', %s, %s" % (key_keyword,key_value,key_comment)
	S.write("""%s
	call FTGTHD(trim(template),card,hdutype,status)
	if (status .ne. 0) then 
		write(*,'(A,I0,A)') "ERROR ",status," PARSING:"
		write(*,'(A)') trim(template)
	endif
	call fl_append_string80_list(cards,card)
"""%template_write)

	value_keyword='DSVAL%.3d'%i
	value_comment=fortran_string(' / Descart parameter %d value' % i,str)
	value_value='options%%%s'%name
	if ptype==str: value_value='trim(%s)'%value_value
	template_write = "	write(template,'(A,%s,A)') '%s  ', %s, %s" % (FORTRAN_CODES[ptype],value_keyword,value_value,value_comment)
	S.write("""%s
	call FTGTHD(template,card,hdutype,status)
	if (status .ne. 0) then 
	write(*,'(A,I0,A)') "ERROR ",status," PARSING:"
		write(*,'(A)') trim(template)
	endif
	call fl_append_string80_list(cards,card)
"""%template_write)


#def fits_card(name,ptype,help,i,S):
#	key_keyword=fortran_string("DSKEY%.3d"%i,str)
#	value_keyword=fortran_string("DSVAL%.3d"%i,str)
#	key_value=fortran_string(name,str)
#	value_value="options%%%s"%name
#	code=FITS_CODES[ptype]
#	
#	#Put a keyword describing the name of the descart parameter:
#	S.write('	call FTPKYS(unit,%s,%s,"Descart parameter %d name",status)\n' % (key_keyword,key_value,i))
#
#	#And then the value
#	if ptype in [str,int,bool]:
#		S.write('	call FTPKY%c(unit,%s,%s,"Descart parameter %d value",status)\n' % (code,value_keyword,value_value,i))
#	elif ptype==float:
#		S.write('	call FTPKYD(unit,%s,%s,3,"Descart parameter %d value",status)\n' % (value_keyword,value_value,i))
#
def generate_parameter_file(options):
	S=StringIO.StringIO()
	S.write("# The default version of this file, %s, is auto-generated by build_options.py.  Edits may be lost when compiling.\n\n"%DEFAULT_PARAMETER_FILENAME)
	for set_name,option_set in options:
		S.write("########################\n")		
		S.write("###### %s\n"%set_name)		
		S.write("########################\n\n")
		for name,ptype,default,help in option_set:
			S.write("#%s\n"%help)
			S.write("%s = %s\n\n"%(name,ini_string(default,ptype)) )
	S.seek(0)
	return S.read()

file(FORTRAN_FILENAME,"w").write(generate_fortran(options))
file(DEFAULT_PARAMETER_FILENAME,"w").write(generate_parameter_file(options))
