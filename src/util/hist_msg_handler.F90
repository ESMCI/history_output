module hist_msg_handler

   !! Class with methods to collect log and error messages for future handling
   !! By default, a <hist_log_messages> object only logs errors
   !! Use the <set_log_level> to change the logging level
   !! The logger stops collecting non-ERROR messages after an error is logged
   !! Use the <new_error>, <add_stack_frame> or <new_message> methods to
   !!   record new messages
   !! Use the <output> method to write the messages to a Fortran unit

   implicit none
   private

   ! Public convenience interfaces
   public :: hist_have_error
   public :: hist_add_message
   public :: hist_add_error
   public :: hist_add_alloc_error

   ! Message collection levels
   !  A lower value is a higher priority
   !  Only collect messages with a level lower or equal to the priority level
   integer, public, parameter :: ERROR   = 0
   integer, public, parameter :: WARNING = 1
   integer, public, parameter :: INFO    = 2
   integer, public, parameter :: VERBOSE = 3
   integer, public, parameter :: DEBUG_HIST   = 4
   integer, public, parameter :: MAX_LEVEL = DEBUG_HIST

   type :: hist_log_entry
      integer                           :: message_level = ERROR
      character(len=:),     allocatable :: log_message
      type(hist_log_entry), pointer     :: next => NULL()
   contains
      final :: finalize_log_entry
   end type hist_log_entry

   type, public :: hist_log_messages
      integer,                private        :: length = 0
      integer,                private        :: error_count = 0
      logical,                private        :: alloc_error = .false.
      integer                                :: message_level = ERROR
      integer,                private        :: max_len = -1 ! no limit
      type(hist_log_entry), private, pointer :: log_messages => NULL()
      type(hist_log_entry), private, pointer :: tail => NULL()
   contains
      procedure :: num_errors      => hist_msgs_num_errors
      procedure :: num_messages    => hist_msgs_num_messages
      procedure :: new_error       => hist_msgs_new_error
      procedure :: new_message     => hist_msgs_new_message
      procedure :: add_alloc_error => hist_msgs_new_alloc_error
      procedure :: add_stack_frame => hist_add_stack_frame
      procedure :: output          => hist_msgs_print_log
      final     :: finalize_hist_log_messages
   end type hist_log_messages

   character(len=8) :: MSG_HEAD(ERROR:DEBUG_HIST) = (/ 'ERROR   ', 'WARNING ',     &
        'INFO    ', 'VERBOSE ', 'DEBUG   ' /)

CONTAINS

   !#######################################################################

   logical function hist_have_error(errors)
      ! Return .true. iff <errors> is present and contains error messages
      type(hist_log_messages), optional, intent(inout) :: errors

      hist_have_error = present(errors)
      if (hist_have_error) then
         hist_have_error = errors%num_errors() > 0
      end if
   end function hist_have_error

   !#######################################################################

   subroutine hist_add_message(subname, msg_level, msgstr1, msgint1, msgstr2, &
        msgint2, msgstr3, logger)
      ! Record a new error message if <errors> is present
      character(len=*),                  intent(in)    :: subname
      integer,                           intent(in)    :: msg_level
      character(len=*),                  intent(in)    :: msgstr1
      integer,                 optional, intent(in)    :: msgint1
      character(len=*),        optional, intent(in)    :: msgstr2
      integer,                 optional, intent(in)    :: msgint2
      character(len=*),        optional, intent(in)    :: msgstr3
      type(hist_log_messages), optional, intent(inout) :: logger

      if (present(logger)) then
         call logger%new_message(msg_level, msgstr1, msgint1=msgint1,         &
              msgstr2=msgstr2, msgint2=msgint2, msgstr3=msgstr3,              &
              subname=subname)
      end if
   end subroutine hist_add_message

   !#######################################################################

   subroutine hist_add_error(subname, errstr1, errint1, errstr2, errint2,     &
        errstr3, errors)
      ! Record a new error message if <errors> is present
      character(len=*),                  intent(in)    :: subname
      character(len=*),                  intent(in)    :: errstr1
      integer,                 optional, intent(in)    :: errint1
      character(len=*),        optional, intent(in)    :: errstr2
      integer,                 optional, intent(in)    :: errint2
      character(len=*),        optional, intent(in)    :: errstr3
      type(hist_log_messages), optional, intent(inout) :: errors

      if (present(errors)) then
         call errors%new_error(errstr1, errint1=errint1, errstr2=errstr2,     &
              errint2=errint2, errstr3=errstr3, subname=subname)
      end if
   end subroutine hist_add_error

   !#######################################################################

   subroutine hist_add_alloc_error(fieldname, filename, line, subname, errors)
      ! Dummy Arguments
      character(len=*),                  intent(in)    :: fieldname
      character(len=*),                  intent(in)    :: filename
      integer,                           intent(in)    :: line
      character(len=*),        optional, intent(in)    :: subname
      type(hist_log_messages), optional, intent(inout) :: errors

      if (present(errors)) then
         call errors%add_alloc_error(fieldname, filename, line, subname=subname)
      end if
   end subroutine hist_add_alloc_error

   !#######################################################################

   subroutine finalize_log_entry(this)
      ! Dummy Argument
      type(hist_log_entry) :: this
      ! Local argument
      type(hist_log_entry), pointer :: next

      if (allocated(this%log_message)) then
         deallocate(this%log_message)
      end if

      next => this%next
      nullify(this%next)
      if (associated(next)) then
         deallocate(next)
      end if
   end subroutine finalize_log_entry

   !#######################################################################

   integer function hist_msgs_num_messages(this)
      ! Dummy Argument
      class(hist_log_messages) :: this

      hist_msgs_num_messages = this%length
   end function hist_msgs_num_messages

   !#######################################################################

   integer function hist_msgs_num_errors(this)
      ! Dummy Argument
      class(hist_log_messages) :: this

      hist_msgs_num_errors = this%error_count
   end function hist_msgs_num_errors

   !#######################################################################

   subroutine hist_msgs_new_error(this, errstr1,                              &
        errint1, errstr2, errint2, errstr3, subname)
      ! Dummy Arguments
      class(hist_log_messages)               :: this
      character(len=*),           intent(in) :: errstr1
      integer,          optional, intent(in) :: errint1
      character(len=*), optional, intent(in) :: errstr2
      integer,          optional, intent(in) :: errint2
      character(len=*), optional, intent(in) :: errstr3
      character(len=*), optional, intent(in) :: subname

      call this%new_message(ERROR, errstr1, msgint1=errint1, msgstr2=errstr2, &
           msgint2=errint2, msgstr3=errstr3, subname=subname)
      this%error_count = this%error_count + 1

   end subroutine hist_msgs_new_error

   !#######################################################################

   subroutine hist_msgs_new_message(this, msg_level, msgstr1,                 &
        msgint1, msgstr2, msgint2, msgstr3, subname)
      ! Dummy Arguments
      class(hist_log_messages)               :: this
      integer,                    intent(in) :: msg_level
      character(len=*),           intent(in) :: msgstr1
      integer,          optional, intent(in) :: msgint1
      character(len=*), optional, intent(in) :: msgstr2
      integer,          optional, intent(in) :: msgint2
      character(len=*), optional, intent(in) :: msgstr3
      character(len=*), optional, intent(in) :: subname

      ! Local variables
      integer,              parameter   :: num_strs = 3
      integer                           :: aerr
      integer                           :: msg_lvl
      integer                           :: msg_len
      type(hist_log_entry), pointer     :: new_log_entry
      character(len=256)                :: msg_strs(num_strs)
      character(len=16)                 :: int_strs(num_strs - 1)

      if ((msg_level < ERROR) .or. (msg_level > MAX_LEVEL)) then
         if (present(subname)) then
            call hist_add_message(subname, WARNING, 'msg_level, ',            &
                 msgint1=msg_level, msgstr2=' out of range, using ERROR',     &
                 logger=this)
         else
            call hist_add_message('hist_msgs_new_message', WARNING,           &
                 'msg_level, ', msgint1=msg_level,                            &
                 msgstr2=' out of range, using ERROR', logger=this)
         end if
         msg_lvl = ERROR
      else
         msg_lvl = msg_level
      end if
      ! Only output if message level is appropriate:
      !   No errors: check message level
      !   Errors: only log errors
      if (((this%num_errors() == 0) .and.                                     &
           (msg_lvl <= this%message_level)) .or.                              &
           (msg_lvl == ERROR)) then
         ! First, calculate the message length
         msg_strs(:) = ''
         int_strs(:) = ''
         msg_len = len_trim(msgstr1)
         if (msg_lvl == ERROR) then
            msg_len = msg_len + len_trim(MSG_HEAD(msg_lvl))
         end if
         if (present(subname)) then
            if (msg_lvl == ERROR) then
               write(msg_strs(1), '(3a)') " in ", trim(subname), ": "
            else
               write(msg_strs(1), '(3a)') "In ", trim(subname), ": "
            end if
         else
            write(msg_strs(1), '(3a)') ": "
         end if
         msg_len = msg_len + len_trim(msg_strs(1)) + 1 ! Keep space at end
         if (present(msgint1)) then
            write(int_strs(1), '(i0)') msgint1
            msg_len = msg_len + len_trim(int_strs(1))
         end if
         if (present(msgint2)) then
            write(int_strs(2), '(i0)') msgint2
            msg_len = msg_len + len_trim(int_strs(2))
         end if
         if (present(msgstr2)) then
            msg_strs(2) = trim(msgstr2)
            msg_len = msg_len + len_trim(msgstr2)
         end if
         if (present(msgstr3)) then
            msg_strs(3) = trim(msgstr3)
            msg_len = msg_len + len_trim(msgstr3)
         end if
         allocate(new_log_entry, stat=aerr)
         ! We can't really report an error here so just trap and record
         if (aerr == 0) then
            allocate(character(len=msg_len) :: new_log_entry%log_message,     &
                 stat=aerr)
            if (aerr == 0) then
               if (msg_lvl == ERROR) then
                  write(new_log_entry%log_message, '(8a)')                    &
                       trim(MSG_HEAD(msg_lvl)), trim(msg_strs(1)), ' ',       &
                       trim(msgstr1), trim(int_strs(1)), trim(msg_strs(2)),   &
                       trim(int_strs(2)), trim(msg_strs(3))
               else
                  write(new_log_entry%log_message, '(7a)') trim(msg_strs(1)), &
                       ' ', trim(msgstr1), trim(int_strs(1)),                 &
                       trim(msg_strs(2)), trim(int_strs(2)), trim(msg_strs(3))
               end if
            else
               ! Hail Mary attempt
               new_log_entry%log_message = 'Message allocation error'
               this%alloc_error = .true.
            end if
            if (.not. associated(this%log_messages)) then
               this%log_messages => new_log_entry
            end if
            if (associated(this%tail)) then
               this%tail%next => new_log_entry
            end if
            this%tail => new_log_entry
            this%length = this%length + 1
         else
            this%alloc_error = .true.
         end if
      end if
   end subroutine hist_msgs_new_message

   !#######################################################################

   subroutine hist_add_stack_frame(this, msg_level, filename, line, subname)
      ! Dummy Arguments
      class(hist_log_messages)               :: this
      integer,                    intent(in) :: msg_level
      character(len=*),           intent(in) :: filename
      integer,                    intent(in) :: line
      character(len=*), optional, intent(in) :: subname

      if (present(subname)) then
         call this%new_message(msg_level,                                     &
              'called from '//trim(subname)//' @'//trim(filename)//':',       &
              msgint1=line)
      else
         call this%new_message(msg_level, 'called from '//trim(filename)//':', &
              msgint1=line)
      end if

   end subroutine hist_add_stack_frame

   !#######################################################################

   subroutine hist_msgs_new_alloc_error(this, fieldname, filename, line,      &
        subname)
      ! Dummy Arguments
      class(hist_log_messages)               :: this
      character(len=*),           intent(in) :: fieldname
      character(len=*),           intent(in) :: filename
      integer,                    intent(in) :: line
      character(len=*), optional, intent(in) :: subname

      if (present(subname)) then
         call this%new_message(ERROR,                                         &
              "Unable to allocate '"//trim(fieldname)//"'",                   &
              msgstr2=' at '//trim(subname)//' @'//trim(filename)//':',       &
              msgint2=line)
      else
         call this%new_message(ERROR,                                         &
              "Unable to allocate '"//trim(fieldname)//"'",                   &
              msgstr2=' at '//trim(filename)//':', msgint2=line)
      end if
   end subroutine hist_msgs_new_alloc_error

   !#######################################################################

   subroutine hist_msgs_print_log(this, unit, header)
      ! Dummy Arguments
      class(hist_log_messages)                   :: this
      integer,                        intent(in) :: unit
      character(len=*),     optional, intent(in) :: header

      type(hist_log_entry), pointer              :: msg_ptr

      msg_ptr => this%log_messages
      if (present(header) .and. (this%num_errors() > 0)) then
         write(unit, '(a)') trim(header)
      end if
      if (this%alloc_error) then
         write(unit, '(a)') 'ERROR: Allocation error in error-handling system'
      end if
      do
         if (associated(msg_ptr)) then
            write(unit, '(a)') trim(msg_ptr%log_message)
            msg_ptr => msg_ptr%next
         else
            exit
         end if
      end do
   end subroutine hist_msgs_print_log

   !#######################################################################

   subroutine finalize_hist_log_messages(this)
      ! Dummy Argument
      type(hist_log_messages) :: this

      if (associated(this%log_messages)) then
         deallocate(this%log_messages)
      end if
      nullify(this%log_messages)
      nullify(this%tail)
      this%length = 0
   end subroutine finalize_hist_log_messages

end module hist_msg_handler
