!! ZM: initital try of land use transitions implementation, 
!!      workflow: in the memory.f90 define the variables, 
!!                in the config.f90 define the namelist
!!                in the init.f90 read in data as a time loop or just define the name of the read in data, and read in the data in the time loop
!!                in the interface.f90 define the connection with the ICON_LAND infrastructure, such as inout and output the the process you are implementing


!! inline comments of Stefanie is in branch feature/landuse_transitions
!! my own developing branch is a outbranch of land4icon_mpim: feature/Landuse_transitions

!! check out the toy model 
MODULE mo_jsb_transitions

    USE mo_util,              ONLY: one_of
    USE mo_jsb_time,          ONLY: is_newday, is_newyear
    USE mo_jsb_process_class, ONLY: Get_process_id
    USE mo_jsb_lcc_class,     ONLY: t_jsb_lcc_proc
    USE mo_jsb_lcc,           ONLY: init_lcc_reloc, start_lcc_reloc, end_lcc_reloc, transfer_active_to_passive_onChunk

! read in data
    landuse_classes = (/TRIM("crop"),TRIM("natl"),TRIM("past")/) ! should be update directly when read in the data
    input_transition_file = alcc_transition_filename_prefix_ + is_newyear

   


    


END MODULE mo_jsb_transitions

