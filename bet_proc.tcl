#

#   bet_proc.tcl - GUI proc for BET - Brain Extraction Tool
#
#   Stephen Smith, FMRIB Image Analysis Group
#
#   Copyright (C) 1999-2000 University of Oxford
#
#   Part of FSL - FMRIB's Software Library
#   WWW:      http://www.fmrib.ox.ac.uk/fsl
#   Email:    fsl@fmrib.ox.ac.uk
#   
#   Developed at FMRIB (Oxford Centre for Functional Magnetic Resonance
#   Imaging of the Brain), Department of Clinical Neurology, Oxford
#   University, Oxford, UK
#   
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or (at
#   your option) any later version.
#   
#   This program is distributed in the hope that it will be useful, but in
#   order that the University as a charitable foundation protects its
#   assets for the benefit of its educational and research purposes, the
#   University makes clear that no condition is made or to be implied, nor
#   is any warranty given or to be implied, as to the accuracy of FSL, or
#   that it will be suitable for any particular purpose or for use under
#   any specific conditions.  Furthermore, the University disclaims all
#   responsibility for the use which is made of FSL.  See the GNU General
#   Public License for more details.
#   
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
#   USA

proc bet_proc { In Out segment_yn overlay_yn mask_yn threshold_yn xtopol_yn cost_yn skull_yn combine_skull_yn fraction gradient } {
    #{{{ setup for running bet 

global PXHOME USER FSLDIR PROCID INMEDX

if { $INMEDX } {
    MxGetCurrentFolder Folder
    MxGetCurrentPage In

    if { [ what_type 0 ] != 3 } {
	puts "Selected page isn't a volume!"
	return 1
    }

    MxGetImageProperties $In InProps
    set input_name [keylget InProps Name]
    
    set InF [ exec ${FSLDIR}/bin/tmpnam /tmp/bet ]

    FSLSaveAs $In AVW ${InF}.hdr true

    set OutF ${InF}_brain
} else {
    set InF  [ file rootname $In ]
    set OutF [ file rootname $Out ]
}

#}}}
    #{{{ run command

set thecommand "exec $FSLDIR/bin/bet $InF $OutF -f $fraction -g $gradient"

if { ! $segment_yn } {
    set thecommand "${thecommand} -n"
}

if { $overlay_yn } {
    set thecommand "${thecommand} -o"
}

if { $mask_yn } {
    set thecommand "${thecommand} -m"
}

if { $threshold_yn } {
    set thecommand "${thecommand} -t"
}

if { $xtopol_yn } {
    set thecommand "${thecommand} -x"
}

if { $cost_yn } {
    set thecommand "${thecommand} -c"
}

if { $skull_yn } {
    set thecommand "${thecommand} -s"
}

if { $INMEDX } {
    ScriptUpdate "Running $thecommand"
}

puts $thecommand

catch { eval $thecommand } ErrMsg

puts "$ErrMsg\nFinished"

if { $INMEDX } {
    CancelScriptUpdate
}

#}}}
    #{{{ read outputs into MEDx

if { $INMEDX } {

    if { $mask_yn } {
	MxOpenImage $Folder ${OutF}_mask.img Mask
	MxGetImageProperties $Mask MaskProps
	keylset MaskProps Name "Brain mask from ${input_name}"
	MxSetImageProperties $Mask $MaskProps
    }

    if { $cost_yn } {
	MxOpenImage $Folder ${OutF}_cost.img Cost
	MxGetImageProperties $Cost CostProps
	keylset CostProps Name "Brain cost from ${input_name}"
	MxSetImageProperties $Cost $CostProps
    }

    if { $skull_yn } {
	MxOpenImage $Folder ${OutF}_skull.img Skull
	MxGetImageProperties $Skull SkullProps
	keylset SkullProps Name "Skull from ${input_name}"
	MxSetImageProperties $Skull $SkullProps
    }

    if { $combine_skull_yn } {
	MxQuietStatistics $In { Minimum Maximum } false Stats
        set origmin [ keylget Stats Minimum ]
        set origmax [ keylget Stats Maximum ]
        render_proc $Folder 0 0 0 1 "$In $origmin $origmax $Skull 50 150"
    }

    if { $xtopol_yn } {
	puts "To run xtopol use\n~flitney/src/XTopol/xtopol -coo ${OutF}.coo -dat ${OutF}.dat"
    }

    if { $overlay_yn } {
	MxOpenImage $Folder ${OutF}_overlay.img Overlay
	MxGetImageProperties $Overlay OverlayProps
	keylset OverlayProps Name "Brain surface from ${input_name}"
	MxSetImageProperties $Overlay $OverlayProps
    }

    if { $segment_yn } {
	MxOpenImage $Folder ${OutF}.img Segment
	MxGetImageProperties $Segment SegmentProps
	keylset SegmentProps Name "Brain from ${input_name}"
	MxSetImageProperties $Segment $SegmentProps
    }

    exec sh -c "rm -f ${InF}* ${OutF}*"
}

#}}}
}
