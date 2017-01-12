#!/usr/bin/wish


###############################
# Program to view QUIET scans # 
###############################

package require Itcl 
package require Itk
package require Iwidgets
namespace import itcl::*

######################
#  Global variables  #
######################

set accept_file "NA"
set datasets    [list]
set num_segments 0

set num_modules  19
set num_diodes    4

######################
#  Utility functions #
######################

proc get_run_id {filename} {
    set elements [split $filename "_."]
    set keyword [lsearch -all -inline -regexp $elements {^seg}]
    set i [expr [lsearch $elements $keyword]-1]
    return [lindex $elements $i]
}

proc get_segment_id {filename} {
    if {[regexp {seg0*([0-9]+)} $filename s t]} {
	return [expr $t]
    } else {
	return -1
    }
}

proc get_module_id {filename} { 
    if {[regexp {mod0*([0-9]+)} $filename s t]} {
	return [expr $t]
    } else {
	return -1
    }
}

proc get_diode_id {filename} {
    if {[regexp {di0*([0-9]+)} $filename s t]} {
	return [expr $t]
    } else {
	return -1
    }
}


proc save_and_exit {} {
    global output_file
    global datasets

    output_database

    toplevel .top 
    label .top.confirmation -text "Are you sure you want to exit?"
    button .top.exit -text "Exit" -command {exit}
    button .top.cancel -text "Cancel" -command { destroy .top }
    
    grid .top.confirmation -in .top -column 1 -row 1 -columnspan 2
    grid .top.exit   -in .top -column 1 -row 2 
    grid .top.cancel -in .top -column 2 -row 2 
}


proc lremove {listVariable value} {
    upvar 1 $listVariable var
    set idx [lsearch -exact $var $value]
    set var [lreplace $var $idx $idx]
}





######################
#  Class definitions #
######################
class Dataset {
    variable parent ""
    variable children ""

    private variable filename
    private variable run
    private variable segment
    private variable module
    private variable diode
    private variable status
    private variable polarization
    private variable current_child
    private variable num_accepted_diodes
    private variable total_num_diodes
    private variable num_modules 
    private variable num_diodes  
    private variable first_module
    private variable num_diodes_per_segment

    constructor {f r s m d p} {
	set filename      $f
	set run           $r
	set segment       $s
	set module        $m
	set diode         $d
	set polarization  $p
	set status        "Accepted"
	set current_child 0
	if {$p == "T"} {
	    set num_modules   2
	    set num_diodes    2
	    set first_module 17
	} else {
	    set num_modules  17
	    set num_diodes    4
	    set first_module  0
	}
	set num_diodes_per_segment [expr $num_modules*$num_diodes]
    }
    destructor {
	clear
    }

    method add {obj} {
	$obj parent $this
	lappend children $obj
    }
    method clear {} {
	if {$children != ""} {
	    eval delete object $children
	}
	set $children ""
    }
    method parent {pobj} {
	set parent $pobj
    }
    method contents {} {
	return $children
    }
    method get_filename {} {
	return $filename
    }
    method get_status {} {
	return $status
    }
    method set_status {new_status} {
	set status $new_status
    }
    method accept {} {
	if {$status == "Rejected"} {
	    set net_change [$this accept_down]
	    $parent accept_up $net_change
	}
    }
    method accept_up { net_change } {
	set status "Accepted"
	set num_accepted_diodes [expr $num_accepted_diodes+$net_change]
	if {$run != -1} {
	    $parent accept_up $net_change
	}
    }
    method accept_down {} {
	if {$status == "Accepted"} {
	    return 0
	}
	
	set status "Accepted"
	set net_change 0
	if {[llength $children] == 0} {
	    set net_change $total_num_diodes
	} else {
	    foreach child $children {
		set sub_change [$child accept_down]
		set net_change [expr $net_change+$sub_change]
	    }
	}
	set num_accepted_diodes [expr $num_accepted_diodes + $net_change]
	return $net_change
    }
    method reject {} {
	if {$status == "Accepted"} {
	    set net_change [$this reject_down]
	    $parent reject_up $net_change
	}
    }
    method reject_up { net_change } {
	set num_accepted_diodes [expr $num_accepted_diodes-$net_change]
	set new_state "Rejected"
	foreach child $children {
	    set child_state [$child get_status]
	    if {$child_state == "Accepted"} {
		set new_state "Accepted"
	    }
	}
	set status $new_state
	if {$run != -1} {
	    $parent reject_up $net_change
	}
    }
    method reject_down {} {
	if {$status == "Rejected"} {
	    return 0
	}
	
	set status "Rejected"
	set net_change 0
	if {[llength $children] == 0} {
	    set net_change $total_num_diodes
	} else {
	    foreach child $children {
		set sub_change [$child reject_down]
		set net_change [expr $net_change+$sub_change]
	    }
	}
	set num_accepted_diodes [expr $num_accepted_diodes - $net_change]
	return $net_change
    }
    method get_run {} {
	return $run
    }
    method get_segment {} {
	return $segment
    }
    method get_module {} {
	return $module
    }
    method get_diode {} {
	return $diode
    }
    method get_pol_state {} {
	return $polarization
    }
    method get_first_child {} {
	if {[llength $children] > 0} {
	    return [lindex $children 0]
	} else {
	    return $this
	}
    }
    method get_current_child {} {
	return [lindex $children $current_child]
    }
    method get_current_child_index {} {
	return $current_child
    }
    method incr_current_child {} {
	if {$current_child < [expr [llength $children]-1]} {
	    set current_child [expr $current_child+1]
	}
    }
    method decr_current_child {} {
	if {$current_child > 0} {
	    set current_child [expr $current_child-1]
	}
    }
    method get_next_sibling {} {
	set p $parent
	$p incr_current_child
	return [$p get_current_child]
    }
    method get_previous_sibling {} {
	set p $parent
	$p decr_current_child
	return [$p get_current_child]
    }
    method get_parent {} {
	if {$module != -1} {
	    return $parent
	} else {
	    return $this
	}
    }
    method get_parent_no_test {} {
	return $parent
    }
    method initialize_total_diode_count {} {
	if {$diode != -1} {
	    set total_num_diodes 1
	} elseif {$module != -1} {
	    set total_num_diodes $num_diodes
	} elseif {$segment != -1} {
	    set total_num_diodes [expr $num_modules*$num_diodes]
	} else {
	    set total_num_diodes [expr $num_modules*$num_diodes*[llength $children]]
	}
	foreach child $children {
	    $child initialize_total_diode_count
	}
    }
    method initialize_accepted_diode_count {} {
	if {[llength $children] == 0} {
	    if {$status == "Accepted"} {
		set num_accepted_diodes $total_num_diodes
	    } else {
		set num_accepted_diodes 0
	    }
	} else {
	    set num_accepted_diodes 0
	    foreach child $children {
		$child initialize_accepted_diode_count
		set num_accepted_diodes [expr $num_accepted_diodes + [$child get_num_accepted_diodes]]
	    }
	}
    }
    method get_total_num_diodes {} {
	return $total_num_diodes
    }
    method get_num_accepted_diodes {} {
	return $num_accepted_diodes
    }
    method get_accepted_ratio {} {
	return [expr $num_accepted_diodes / (1.*$total_num_diodes)]
    }
    method get_num_accepted_CES {} {
	return [expr $num_accepted_diodes / (1.*$num_diodes_per_segment)]
    }
    method get_node {r s m d} {
	if {$diode != -1} {
	    if {$r == $run && $s == $segment && $m == $module && $d == $diode} {
		return $this
	    } else {
		return -1
	    }
	} elseif {$m == -1 && $module != -1} {
	    if {$r == $run && $s == $segment && $m == $module} {
		return $this
	    } else {
		return -1
	    }
	} else {
	    foreach child $children {
		if {    [expr $r == [$child get_run]     || [expr $run     == -1]] && \
			[expr $s == [$child get_segment] || [expr $segment == -1]] && \
			[expr $m == [$child get_module]  || [expr $module  == -1]]} {
		    set candidate [$child get_node $r $s $m $d]
		    if {$candidate != -1} {
			return $candidate
		    }
		}
	    }   
	    return -1
	}
    }
    method get_num_segments {} {
	return [llength $children]
    }
    method get_segment_obj {r s} {
	foreach child $children {
	    if {$r == [$child get_run] && $s == [$child get_segment]} {
		return $child
	    }
	}
	return -1
    }
    method get_module_obj {r s m} {
	foreach child $children {
	    if {$r == [$child get_run] && $s == [$child get_segment] && $m == [$child get_module]} {
		return $child
	    }
	}
	return -1
    }
    method get_diode_obj {r s m d} {
	foreach child $children {
	    if {$r == [$child get_run] && $s == [$child get_segment] && $m == \
		    [$child get_module] && $d == [$child get_diode]} {
		return $child
	    }
	}
	return -1
    }
    method output_children {} {
	puts "Current object $filename"
	foreach child $children {
	    puts "    Child = [eval $child get_filename]"
	}
    }
    method output_accept_file {out} {
	foreach child $children {
	    # Create a rejected list with everything on by default
	    set reject_list [list]
	    for {set j 0} {$j < $num_modules} {incr j} {
		set m [expr $j+$first_module]
		for {set i 0} {$i < 4} {incr i} {
		    lappend reject_list "$m $i"
		}
	    }	    

	    set reject_list [$child get_accepted_diodes $reject_list]
	    set num_reject [expr $total_num_diodes-$num_accepted_diodes]
	    puts -nonewline $out "[$child get_run] [$child get_segment] [llength $reject_list]"
	    foreach reject $reject_list {
		puts -nonewline $out " $reject"
	    }
	    puts $out ""
	}
    }
    method get_accepted_diodes {reject_list} {
	if {[llength $children] == 0 && $status == "Accepted"} {
	    if {$diode != -1} {
		lremove reject_list "$module $diode" 
		return $reject_list
	    } elseif {$module != -1} {
		for {set i 0} {$i < $num_diodes} {incr i} {
		    lremove reject_list "$module $i" 
		}
		return $reject_list
	    } else {
		for {set j 0} {$j < $num_modules} {incr j} {
		    set m [expr $j+$first_module]
		    for {set i 0} {$i < $num_diodes} {incr i} {
			lremove reject_list "$m $i" 
		    }
		}
		return $reject_list
	    }
	} else {
	    foreach child $children {
		set reject_list [$child get_accepted_diodes $reject_list]
	    }
	    return $reject_list
	}
    }
}



##################
# Initialization # 
##################

# Read command line arguments
set num_datasets 0
for {set i 0} {$i < $argc} {incr i} {
    set arg [lindex $argv $i]
    if {$arg == "-accept" || $arg == "-a"}  {
	set accept_file [lindex $argv [expr $i+1]]
	incr i
    } else {
	lappend datasets [Dataset #auto $arg -1 -1 -1 -1 "T"]
	incr num_datasets
	lappend datasets [Dataset #auto $arg -1 -1 -1 -1 "P"]
	incr num_datasets
    }
}
if {$num_datasets == 0} {
    lappend datasets [Dataset #auto "." -1 -1 -1 -1 "T"]
    incr num_datasets
    lappend datasets [Dataset #auto "." -1 -1 -1 -1 "P"]
    incr num_datasets
}

# Initialize file lists
foreach s $datasets {
    # Check that directory exist
    set directory [eval $s get_filename]
    if {[file exists $directory]} {
	set pol [eval $s get_pol_state]
	set dirlist [glob -nocomplain "$directory/*.gif"]
	set filelist [list]
	foreach file $dirlist {
	    if { [regexp "seg\[0-9\]*_${pol}(_map\|_rhs\|_bin).gif"  $file]}  {
		lappend filelist $file
	    }	    
	}
	if {[llength $filelist] == 0} {
	    lremove datasets $s
	} else {
	    set filelist [lsort -dictionary $filelist]
	    set numfiles [llength $filelist]
	    puts "Total number of files = $numfiles"    

	    # Add segment files
	    foreach file $filelist {
		if { [regexp "seg\[0-9\]*_${pol}(_map\|_rhs\|_bin).gif"  $file]}  {
		    $s add [Dataset #auto $file [get_run_id $file] [get_segment_id $file] -1 -1 $pol]
		}
	    }

	    # Add module files
	    foreach file $filelist {
		if { [regexp "mod\[0-9\]*_${pol}(_map\|_rhs\|_bin).gif"  $file]}  {
		    set seg [$s get_segment_obj [get_run_id $file] [get_segment_id $file]]
		    if {$seg != -1} {
			$seg add [Dataset #auto $file [get_run_id $file] [get_segment_id $file] \
				      [get_module_id $file] -1 $pol]
		    }
		}
	    }

	    # Add diode files
	    foreach file $filelist {
		if { [regexp "di\[0-9\]*_${pol}(_map\|_rhs\|_bin).gif"  $file]}  {
		    set seg [$s get_segment_obj [get_run_id $file] [get_segment_id $file]]
		    set mod [$seg get_module_obj [get_run_id $file] [get_segment_id $file] [get_module_id $file]]
		    if {$mod != -1} {
			$mod add [Dataset #auto $file [get_run_id $file] [get_segment_id $file] \
				      [get_module_id $file] [get_diode_id $file] $pol]
		    }
		}
	    }
	    
	}
    } else {
	lremove datasets $s
    }
}

# Initialize diode counts
foreach s $datasets {
    $s initialize_total_diode_count
    $s initialize_accepted_diode_count
}

#  Process accept file
if {$accept_file != "NA" && [file exists $accept_file]} {
    set fp [open "$accept_file" r]
    set data [read $fp]
    close $fp

    set data [split $data "\n"]
    foreach line $data {
	set elements [split $line]
	if {$elements != ""} {
	    set run        [lindex $elements 0]
	    set segment    [lindex $elements 1]
	    set num_reject [lindex $elements 2]
	    for {set i 0} {$i < [expr 2*$num_reject]} {incr i 2} {
		set module [lindex $elements [expr 3+$i]]
		set diode  [lindex $elements [expr 4+$i]]

		# Reject from all data sets
		foreach s $datasets {
		    set node [$s get_node $run $segment $module $diode]
		    if {$node != -1} {
			$node reject
#			$node set_status "Rejected"
		    }
		}
	    }
	}	
    }
}


####################
# Global variables # 
####################

set current_dataset       [lindex $datasets 0]
set current_obj           [$current_dataset get_first_child]
set current_filenum       1
set current_file          ""
set current_status        "Accepted"
set current_num_segments  0
set num_accept_CES        0
set accept_ratio          0
set L2_directory          "/xanadu/project/quiet"
set random_number         [expr rand()]
set todfile               "/tmp/tod_${random_number}.dat"
set psfile                "/tmp/ps_${random_number}.dat"
set viewer                "gnuplot"


################
# GUI elements # 
################

proc update_status {} {
    global current_obj
    global current_filenum
    global current_file
    global current_status
    global current_num_segments
    global num_accept_CES
    global accept_ratio
    global current_dataset
    set current_file    [$current_obj get_filename]
    set current_status  [$current_obj get_status]
    if {[$current_obj get_module] == -1} {
	set p [$current_obj get_parent_no_test]
	set current_num_segments [$p get_num_segments]
	set current_filenum      [expr [$p get_current_child_index]+1]
    }
    set num_accept_CES [$current_dataset get_num_accepted_CES]
    set accept_ratio   [$current_dataset get_accepted_ratio]
    image create photo imgobj -file [$current_obj get_filename]
}

# Navigation buttons
proc push_next {} {
    global current_obj
    set current_obj [$current_obj get_next_sibling]
    update_status
}

proc push_previous {} {
    global current_obj
    set current_obj [$current_obj get_previous_sibling]
    update_status
}

proc jump_up {} {
    global current_obj
    set current_obj [$current_obj get_parent]
    update_status
}

proc jump_down {} {
    global current_obj
    set current_obj [$current_obj get_first_child]
    update_status
}

proc reject_element {} {
    global current_obj
    $current_obj reject
    update_status
}

proc accept_element {} {
    global current_obj
    $current_obj accept 
    update_status
}

proc output_database {} {
    global datasets

    foreach s $datasets {
	set filename "accept_$s.txt"
	puts "Outputting accept file = $filename"

	# Write data base
	set out [open $filename w]
	$s output_accept_file $out
	close $out
    }
}

proc basename {path} {
    regexp {.*/(.*)} $path foo filename
    return $filename
}

proc plot_tod { function } {
    global current_obj
    global L2_directory
    global bindir
    global todfile
    global viewer

    set module [$current_obj get_module]
    set diode [$current_obj get_diode]
    if {$diode != -1} {

	set filename [$current_obj get_filename]
	set filename [basename $filename]
	if {[regexp {^([^_]+)_(.*?)(_sub[0-9]+)?_[0-9]*(_seg[0-9]+(_mod[0-9]+(_di[0-9]+)?)?)?_[TP](_map\|_rhs\|_bin).gif} \
		 $filename a patch_type patch_name]} {
	    set l2file "${L2_directory}/${patch_type}/${patch_name}/${patch_type}_${patch_name}_[$current_obj get_run]_seg[format %03d [$current_obj get_segment]].unf"
	    if {[file exists $l2file]} {
		puts "Getting TOD for $l2file"
		exec l2txt $l2file demod\[$module\]\[$diode\] > $todfile
		if {$function == "ps"} {
		    exec powspec $todfile $todfile
		}
		if {$viewer == "gnuplot"} {
		    if {$function == "ps"} {
			exec echo "set logscale xy; plot '$todfile' w l" | gnuplot -persist &
		    } else {
			exec echo "plot '$todfile' w l" | gnuplot -persist &
		    }
		} elseif {$viewer == "xmgrace"} {
		    exec xmgrace $todfile &
		}

	    } else {
		puts "L2 file = $l2file does not exist"
	    }
	}
    }
}




proc push_exit {} {
    toplevel .top 
    label .top.confirmation -text "Are you sure you want to exit?"
    button .top.exit -text "Exit" -command {exit}
    button .top.cancel -text "Cancel" -command { destroy .top }
    
    grid .top.confirmation -in .top -column 1 -row 1 -columnspan 2
    grid .top.exit   -in .top -column 1 -row 2 
    grid .top.cancel -in .top -column 2 -row 2 
}
button .exit_button -text "Exit" -command push_exit

# Canvas
image create photo imgobj -file [$current_obj get_filename]
canvas .image_canvas -relief sunken -background "blue" -height [image height imgobj] -width [image width imgobj]
.image_canvas create image 0 0 -image imgobj -anchor nw


###############
# Set up GUI  # 
###############

update_status



#entry .numfiles_entry -textvariable numfiles -state readonly -width 8

frame .counter_frame 
frame .info_frame 

#grid .curr_text  -in .info_frame -column 1 -row 1 
#grid .numfiles_label -in .counter_frame -column 2 -row 1 

#grid .current_entry  -in .counter_frame -column 1 -row 2 
#grid .numfiles_entry -in .counter_frame -column 2 -row 2 



#label .level_label    -text "Level" 
label .curr_label        -text "Current" 
label .accept_label      -text "Accepted" 
label .numseg_label      -text "Total" 
label .ratio_label       -text "Ratio" 
#label .module_label   -text "Module" 
#label .diode_label    -text "Diode" 
#label .chisq_label    -text "Chi-square" 
label .state_label    -text "State" 

#entry .level_entry    -textvariable current_level    -state readonly -width 8
entry .curr_entry     -textvariable current_filenum      -state readonly -width 8
entry .numseg_entry   -textvariable current_num_segments -state readonly -width 8
entry .accept_entry   -textvariable num_accept_CES       -state readonly -width 8
entry .ratio_entry    -textvariable accept_ratio         -state readonly -width 8
#entry .run_entry      -textvariable current_run      -state readonly -width 8
#entry .segment_entry  -textvariable current_segment  -state readonly -width 8
#entry .module_entry   -textvariable current_module   -state readonly -width 8
#entry .diode_entry    -textvariable current_diode    -state readonly -width 8
#entry .chisq_entry    -textvariable chisq            -state readonly -width 8
entry .state_entry    -textvariable current_status    -state readonly -width 10

#grid .level_label    -in .counter_frame -column 1 -row 1 
grid .curr_label      -in .counter_frame -column 1 -row 1 
grid .numseg_label    -in .counter_frame -column 2 -row 1 
grid .accept_label    -in .counter_frame -column 3 -row 1 
grid .ratio_label     -in .counter_frame -column 4 -row 1 
#grid .run_label      -in .counter_frame -column 4 -row 1 
#grid .segment_label  -in .counter_frame -column 5 -row 1 
#grid .module_label   -in .counter_frame -column 6 -row 1 
#grid .diode_label    -in .counter_frame -column 7 -row 1 
#grid .chisq_label    -in .counter_frame -column 8 -row 1 
grid .state_label    -in .counter_frame -column 5 -row 1 

#grid .level_entry    -in .counter_frame -column 1 -row 2
grid .curr_entry      -in .counter_frame -column 1 -row 2
grid .numseg_entry    -in .counter_frame -column 2 -row 2 
grid .accept_entry    -in .counter_frame -column 3 -row 2
grid .ratio_entry     -in .counter_frame -column 4 -row 2
#grid .run_entry      -in .counter_frame -column 4 -row 2 
#grid .segment_entry  -in .counter_frame -column 5 -row 2 
#grid .module_entry   -in .counter_frame -column 6 -row 2 
#grid .diode_entry    -in .counter_frame -column 7 -row 2 
#grid .chisq_entry    -in .counter_frame -column 8 -row 2 
grid .state_entry    -in .counter_frame -column 5 -row 2 



label .filename_label -text "Filename =" 
entry .filename_entry -textvariable current_file -state readonly -width 60

#label .state_label -text "Current state =" 
#entry .state_entry -textvariable current_state -state readonly

grid .filename_label -in .info_frame -column 1 -row 1 -sticky e 
grid .filename_entry -in .info_frame -column 2 -row 1 -sticky w

#grid .state_label -in .info_frame -column 1 -row 2 -sticky e
#grid .state_entry -in .info_frame -column 2 -row 2 -sticky w

grid .counter_frame -in .info_frame -column 1 -row 3 -columnspan 2

#grid .l2_dir_label -in .info_frame -column 1 -row 5 -sticky e
#grid .l2_dir_entry -in .info_frame -column 2 -row 5 -sticky w

#grid .accept_file_label -in .info_frame -column 1 -row 6 -sticky e
#grid .accept_file_entry -in .info_frame -column 2 -row 6 -sticky w

#grid .reject_file_label -in .info_frame -column 1 -row 7 -sticky e
#grid .reject_file_entry -in .info_frame -column 2 -row 7 -sticky w

#grid .undecided_file_label -in .info_frame -column 1 -row 8 -sticky e
#grid .undecided_file_entry -in .info_frame -column 2 -row 8 -sticky w


###############
# Create menu # 
###############

#Declare that there is a menu
menu .mbar
. config -menu .mbar

#The Main Buttons
.mbar add cascade -label "File" -underline 0 \
      -menu [menu .mbar.file -tearoff 0]
.mbar add cascade -label "Help" -underline 0 \
      -menu [menu .mbar.help -tearoff 0]



proc about_dialog {} {

    toplevel .about
    text .about.txt -width 80 -height 30
    
    .about.txt delete 1.0 end
    .about.txt insert end {
	About Quiet Scan Viewer version 0.2
	----------
	This is a simple viewer for QUIET scans. 
	Written by Hans Kristian Eriksen	
	----------

	Usage:
	    Right arrow = Go to next image
	    Left arrow  = Go to previous image
	    Up arrow    = Go to higher level image
	    Down arrow  = Go to lower level image
	    r           = Reject segment and children
	    a           = Accept segment and children
	    s           = Save current accept file
	    q           = Save and exit
	    t           = Show TOD in external viewer
	    p           = Show power spectrum in external viewer
	    g           = Select gnuplot as external viewer
	    x           = Select xmgrace as external viewer
    }

    button .about.exit -text "Continue" -command {destroy .about}

    grid .about.txt -in .about -row 1 -column 1
    grid .about.exit -in .about -row 2 -column 1
}

## Help ##
set m .mbar.help
$m add command -label "About" -command about_dialog


##############
# Set up GUI # 
##############

bind . <Up>    { jump_up }
bind . <Down>  { jump_down }
bind . <Left>  { push_previous }
bind . <Right> { push_next }
bind . <r>     { reject_element }
bind . <a>     { accept_element }
bind . <s>     { output_database }
bind . <q>     { save_and_exit }
bind . <t>     { plot_tod "tod" }
bind . <p>     { plot_tod "ps" }
bind . <g>     { set viewer "gnuplot"}
bind . <x>     { set viewer "xmgrace"}
bind . <h>     { about_dialog}


# Set up user interface
grid .info_frame   -in . -row 1 -column 1 -columnspan 4
grid .image_canvas -in . -row 2 -column 1 -columnspan 4





