module hist_buffer
   use, intrinsic :: ISO_FORTRAN_ENV, only: REAL64
   use hist_hashable,                 only: hist_hashable_t

   implicit none
   private

   ! Public interfaces
   public :: buffer_factory

   ! Accumulation types -- the integers and array positions below must match
   integer, parameter, public :: hist_accum_lst = 1 ! last sample
   integer, parameter, public :: hist_accum_min = 2 ! minimum sample
   integer, parameter, public :: hist_accum_max = 3 ! maximum sample
   integer, parameter, public :: hist_accum_avg = 4 ! sample average
   integer, parameter, public :: hist_accum_var = 5 ! sample standard deviation

   ! Number of accumulation types defined above
   integer, parameter, public :: hist_num_accum_types = 5

   integer, parameter         :: as_len = 36
   character(len=as_len), parameter, public ::                                &
        accum_strings(hist_num_accum_types) = [                               &
        'last sampled value                  ',                               &
        'minimum of sampled values           ',                               &
        'maximum of sampled values           ',                               &
        'average of sampled values           ',                               &
        'standard deviation of sampled values']

   character(len=3), parameter, public ::                                     &
        accum_abbrev(hist_num_accum_types) =                                  &
        ['lst', 'min', 'max', 'avg', 'var']

   type, abstract, public :: hist_buffer_t
      ! hist_buffer_t is an abstract base class for hist_outfld buffers
      class(hist_hashable_t), pointer              :: field_info => null()
      integer,                             private :: vol = -1 ! For host output
      integer,                             private :: horiz_axis_ind = 0
      integer,                             private :: rank = 0
      integer,                             private :: accum_type = 0
      integer,                allocatable, private :: field_shape(:)
      integer,                allocatable, private :: block_begs(:)
      integer,                allocatable, private :: block_ends(:)
      character(len=:),       allocatable, private :: buff_type
      class(hist_buffer_t),   pointer              :: next => null()
   contains
      procedure                              :: field  => get_field_info
      procedure                              :: volume => get_volume
      procedure                              :: get_shape
      procedure                              :: horiz_axis_index
      procedure                              :: init_buffer
      procedure                              :: accum_string
      procedure                              :: buffer_type
      procedure                              :: check_status
      procedure                              :: has_blocks
      procedure(hist_buff_clear),   deferred :: clear
      procedure(hist_buff_init),    deferred :: initialize
   end type hist_buffer_t

   type, public, extends(hist_buffer_t) :: hist_buff_1d_t
      real(REAL64), allocatable :: data(:)
      real(REAL64), allocatable :: var_buffer(:)
      integer,  allocatable,   private :: num_samples(:)
   contains
      procedure :: clear => buff_1d_clear
      procedure :: accumulate => buff_1d_accum
      procedure :: norm_value => buff_1d_value
      procedure :: initialize => init_buff_1d
   end type hist_buff_1d_t

   type, public, extends(hist_buffer_t) :: hist_buff_2d_t
      real(REAL64), allocatable :: data(:,:)
      real(REAL64), allocatable :: var_buffer(:,:)
      integer, allocatable,    private :: num_samples(:)
   contains
      procedure :: clear => buff_2d_clear
      procedure :: accumulate => buff_2d_accum
      procedure :: norm_value => buff_2d_value
      procedure :: initialize => init_buff_2d
      procedure :: check_fill_value => buff_2d_check_fill
   end type hist_buff_2d_t

   ! Abstract interfaces for hist_buffer_t class
   abstract interface
      subroutine hist_buff_sub_log(this, logger)
         use hist_msg_handler, only: hist_log_messages
         import                   :: hist_buffer_t
         implicit none
         class(hist_buffer_t),              intent(inout) :: this
         type(hist_log_messages), optional, intent(inout) :: logger
      end subroutine hist_buff_sub_log
   end interface

   abstract interface
      subroutine hist_buff_init(this, field_in, volume_in, horiz_axis_in,     &
           accum_type_in, shape_in, block_sizes_in, block_ind_in, logger)
         use hist_msg_handler,             only: hist_log_messages
         use, intrinsic :: ISO_FORTRAN_ENV, only: REAL64
         import                   :: hist_buffer_t
         import                   :: hist_hashable_t
         implicit none
         class(hist_buffer_t),              intent(inout) :: this
         class(hist_hashable_t),  pointer,  intent(in)    :: field_in
         integer,                           intent(in)    :: volume_in
         integer,                           intent(in)    :: horiz_axis_in
         integer,                           intent(in)    :: accum_type_in
         integer,                           intent(in)    :: shape_in(:)
         integer,                 optional, intent(in)    :: block_sizes_in(:)
         integer,                 optional, intent(in)    :: block_ind_in
         type(hist_log_messages), optional, intent(inout) :: logger
      end subroutine hist_buff_init
   end interface

   abstract interface
      subroutine hist_buff_clear(this, logger)
         use hist_msg_handler, only: hist_log_messages
         import                  :: hist_buffer_t
         implicit none
         class(hist_buffer_t),              intent(inout) :: this
         type(hist_log_messages), optional, intent(inout) :: logger
      end subroutine hist_buff_clear
   end interface

contains

   !#######################################################################

   function get_field_info(this) result(field_info)
      class(hist_buffer_t), intent(in) :: this
      class(hist_hashable_t), pointer  :: field_info

      field_info => this%field_info
   end function get_field_info

   !#######################################################################

   pure function get_volume(this) result(volume)
      class(hist_buffer_t), intent(in) :: this
      integer                          :: volume

      volume = this%vol
   end function get_volume

   !#######################################################################

   pure function get_shape(this) result(fshape)
      class(hist_buffer_t), intent(in) :: this
      integer, allocatable :: fshape(:)

      allocate(fshape(size(this%field_shape)))

      fshape = this%field_shape
   end function get_shape

   !#######################################################################

   pure function horiz_axis_index(this) result(axis_index)
      class(hist_buffer_t), intent(in) :: this
      integer                          :: axis_index

      axis_index = this%horiz_axis_ind
   end function horiz_axis_index
   !#######################################################################

   function check_status(this, logger, filename, line) result(is_valid)
      ! Check to see if this buffer is properly initialized
      use hist_msg_handler, only: hist_log_messages, hist_add_error, ERROR

      ! Dummy arguments
      class(hist_buffer_t),              intent(inout) :: this
      type(hist_log_messages), optional, intent(inout) :: logger
      character(len=*),        optional, intent(in)    :: filename
      integer,                 optional, intent(in)    :: line
      logical                                          :: is_valid
      ! Local variable
      character(len=*), parameter :: subname = 'check_status'

      is_valid = .true.
      if ((this%horiz_axis_index() < 1)      .or.                             &
           (.not. allocated(this%field_shape))) then
         is_valid = .false.
         call hist_add_error(subname,                                         &
              "buffer not properly initialized '", errors=logger)
         if (present(filename) .and. present(line) .and. present(logger)) then
            call logger%add_stack_frame(ERROR, filename, line)
         end if
      end if

   end function check_status

   !#######################################################################

   pure function has_blocks(this) result(blocked)
      ! Dummy argument
      class(hist_buffer_t),   intent(in)    :: this
      logical                               :: blocked

      blocked = allocated(this%block_begs) .and. allocated(this%block_ends)

   end function has_blocks

   !#######################################################################

   subroutine init_buffer(this, field_in, volume_in, horiz_axis_in,           &
        accum_type_in, shape_in, block_sizes_in, block_ind_in, logger)
      use hist_msg_handler, only: hist_log_messages, hist_add_alloc_error
      use hist_msg_handler, only: hist_add_error

      ! Dummy arguments
      class(hist_buffer_t),              intent(inout) :: this
      class(hist_hashable_t),  pointer,  intent(in)    :: field_in
      integer,                           intent(in)    :: volume_in
      integer,                           intent(in)    :: horiz_axis_in
      integer,                           intent(in)    :: accum_type_in
      integer,                           intent(in)    :: shape_in(:)
      integer,                 optional, intent(in)    :: block_sizes_in(:)
      integer,                 optional, intent(in)    :: block_ind_in
      type(hist_log_messages), optional, intent(inout) :: logger
      ! Local variables
      integer                     :: astat
      integer                     :: hsize
      character(len=256)          :: errmsg
      character(len=*), parameter :: subname = 'init_buffer'

      ! Sanity check
      this%rank = SIZE(shape_in, 1)
      if ((horiz_axis_in < 1) .or. (horiz_axis_in > this%rank)) then
         call hist_add_error(subname,                                         &
              "horiz_axis_in must be between 1 and '", errint1=this%rank,     &
              errors=logger)
      else
         this%field_info => field_in
         this%vol = volume_in
         this%horiz_axis_ind = horiz_axis_in
         this%accum_type = accum_type_in
         allocate(this%field_shape(size(shape_in, 1)), stat=astat, errmsg=errmsg)
         if (astat == 0) then
            this%field_shape(:) = shape_in(:)
         else
            call hist_add_alloc_error('field_shape, error: '//errmsg, __FILE__, __LINE__ - 4,  &
                 subname=subname, errors=logger)
         end if
      end if
   end subroutine init_buffer

   !#######################################################################

   function accum_string(this) result(ac_str)
      class(hist_buffer_t), intent(in)   :: this
      character(len=as_len), allocatable :: ac_str

      ac_str = accum_strings(this%accum_type)

   end function accum_string

   !#######################################################################

   function buffer_type(this) result(bt_str)
      class(hist_buffer_t), intent(in) :: this
      character(len=:), allocatable    :: bt_str

      if (allocated(this%buff_type)) then
         bt_str = this%buff_type
      else
         bt_str = 'UNKNOWN'
      end if

   end function buffer_type

   !#######################################################################

   subroutine buff_1d_clear(this, logger)
      use hist_msg_handler, only: hist_log_messages, hist_add_alloc_error

      ! Dummy arguments
      class(hist_buff_1d_t),       intent(inout) :: this
      type(hist_log_messages), optional, intent(inout) :: logger
      ! Local variables
      integer                     :: aerr
      character(len=256)          :: errmsg
      character(len=*), parameter :: subname = 'buff_1d_clear'

      if (.not. allocated(this%data)) then
         allocate(this%data(this%field_shape(1)), stat=aerr, errmsg=errmsg)
         if (aerr /= 0) then
            call hist_add_alloc_error('data, error: '//errmsg, __FILE__, __LINE__ - 1,  &
                 subname=subname, errors=logger)
         end if
      end if
      if (.not. allocated(this%var_buffer) .and. this%accum_type == hist_accum_var) then
         allocate(this%var_buffer(this%field_shape(1)), stat=aerr, errmsg=errmsg)
         if (aerr /= 0) then
            call hist_add_alloc_error('var_buffer, error: '//errmsg, __FILE__, __LINE__ - 2, &
                 subname=subname, errors=logger)
         end if
      end if
      if (.not. allocated(this%num_samples)) then
         allocate(this%num_samples(this%field_shape(1)), stat=aerr, errmsg=errmsg)
         if (aerr /= 0) then
            call hist_add_alloc_error('num_samples, error: '//errmsg, __FILE__, __LINE__ - 1, &
                 subname=subname, errors=logger)
         end if
      end if
      this%data = 0.0_REAL64
      this%num_samples = 0
      if (allocated(this%var_buffer)) then
         this%var_buffer = 0.0_REAL64
      end if

   end subroutine buff_1d_clear

   !#######################################################################

   subroutine init_buff_1d(this, field_in, volume_in, horiz_axis_in, &
        accum_type_in, shape_in, block_sizes_in, block_ind_in, logger)
      use hist_msg_handler, only: hist_log_messages

      class(hist_buff_1d_t),             intent(inout) :: this
      class(hist_hashable_t),  pointer,  intent(in)    :: field_in
      integer,                           intent(in)    :: volume_in
      integer,                           intent(in)    :: horiz_axis_in
      integer,                           intent(in)    :: accum_type_in
      integer,                           intent(in)    :: shape_in(:)
      integer,                 optional, intent(in)    :: block_sizes_in(:)
      integer,                 optional, intent(in)    :: block_ind_in
      type(hist_log_messages), optional, intent(inout) :: logger

      call init_buffer(this, field_in, volume_in, horiz_axis_in, accum_type_in, &
              shape_in, block_sizes_in, block_ind_in, logger=logger)
      call this%clear(logger=logger)
      this%buff_type = 'buff_1d'

   end subroutine init_buff_1d

   !#######################################################################

   subroutine buff_1d_accum(this, field, cols_or_block, flag_xyfill, &
         fill_value, cole, logger)
      use hist_msg_handler, only: hist_log_messages, hist_add_error
      ! Dummy arguments
      class(hist_buff_1d_t),       intent(inout) :: this
      real(REAL64),                      intent(in)    :: field(:)
      integer,                           intent(in)    :: cols_or_block
      logical,                           intent(in)    :: flag_xyfill
      real(REAL64),                      intent(in)    :: fill_value
      integer,                 optional, intent(in)    :: cole
      type(hist_log_messages), optional, intent(inout) :: logger
      ! Local variables
      integer            :: col_beg_use
      integer            :: col_end_use
      integer            :: ind1
      real(REAL64)       :: fld_val, tmp
      character(len=256) :: errmsg
      character(len=*), parameter :: subname = 'buff_1d_accum'

      if (this%has_blocks()) then
         col_end_use = this%block_ends(cols_or_block)
      else
         ! Non blocked, <cols_or_block> is the first column index
         col_beg_use = cols_or_block
         if (present(cole)) then
            col_end_use = cole
         else
            col_end_use = col_beg_use +                                       &
                 this%field_shape(this%horiz_axis_ind) - 1
         end if
      end if

      select case (this%accum_type)
      case (hist_accum_lst)
         do ind1 = col_beg_use, col_end_use
            fld_val = field(ind1 - col_beg_use + 1)
            if (flag_xyfill) then
               if (fld_val == fill_value) then
                  ! Set num samples to 0 if this is a fill value
                  ! This will trigger logic in _value routine
                  this%data(ind1) = fill_value
                  this%num_samples(ind1) = 0
               else
                  this%data(ind1) = fld_val
                  this%num_samples(ind1) = 1
               end if
            else
               this%data(ind1) = fld_val
               this%num_samples(ind1) = 1
            end if
         end do
      case (hist_accum_min)
         do ind1 = col_beg_use, col_end_use
            fld_val = field(ind1 - col_beg_use + 1)
            if (flag_xyfill) then
               ! If the buffer (possibly) contains fill values,
               ! only check for minimum if not a fill value
               if (this%num_samples(ind1) == 0) then
                  if (fld_val /= fill_value) then
                     this%data(ind1) = fld_val
                     this%num_samples(ind1) = 1
                  else
                     ! Set to large positive number if this is a fill value
                     ! Also do not change num_samples - if num_samples
                     ! is zero at the end of the accumulation period,
                     ! the value will be overwritten to fill_value
                     this%data(ind1) = HUGE(REAL64)
                  end if
               else if (fld_val < this%data(ind1) .and. fld_val /= fill_value) then
                  this%data(ind1) = fld_val
                  this%num_samples(ind1) = 1
               end if ! No else, we already have the minimum value for this col
            else
               if (this%num_samples(ind1) == 0 .or. fld_val < this%data(ind1)) then
                  this%data(ind1) = fld_val
                  this%num_samples(ind1) = 1
               end if ! No else, we already have the minimum value for this col
            end if
         end do
      case (hist_accum_max)
         do ind1 = col_beg_use, col_end_use
            fld_val = field(ind1 - col_beg_use + 1)
            if (flag_xyfill) then
               ! If the buffer (possibly) contains fill values,
               ! only check for maximum if not a fill value
               if (this%num_samples(ind1) == 0) then
                  if (fld_val /= fill_value) then
                     this%data(ind1) = fld_val
                     this%num_samples(ind1) = 1
                  else
                     ! Set to large negative number if this is a fill value
                     ! Also do not change num_samples - if num_samples
                     ! is zero at the end of the accumulation period,
                     ! the value will be overwritten to fill_value
                     this%data(ind1) = -HUGE(REAL64)
                  end if
               else if (fld_val > this%data(ind1) .and. fld_val /= fill_value) then
                  this%data(ind1) = fld_val
                  this%num_samples(ind1) = 1
               end if ! No else, we already have the maximum value for this col
            else
               if (this%num_samples(ind1) == 0 .or. fld_val > this%data(ind1)) then
                  this%data(ind1) = fld_val
                  this%num_samples(ind1) = 1
               end if ! No else, we already have the maximum value for this col
            end if
         end do
      case (hist_accum_avg)
         do ind1 = col_beg_use, col_end_use
            fld_val = field(ind1 - col_beg_use + 1)
            if (flag_xyfill) then
               ! Only include sample if not the fill value
               if (fld_val /= fill_value) then
                  this%data(ind1) = this%data(ind1) + fld_val
                  this%num_samples(ind1) = this%num_samples(ind1) + 1
               end if
            else
               this%data(ind1) = this%data(ind1) + fld_val
               this%num_samples(ind1) = this%num_samples(ind1) + 1
            end if
         end do
      case (hist_accum_var)
         ! Standard deviation using Welford's algorithm
         ! DOI: 10.1080/00401706.1962.10490022
         do ind1 = col_beg_use, col_end_use
            fld_val = field(ind1 - col_beg_use + 1)
            if (flag_xyfill) then
               if (fld_val /= fill_value) then
                  ! Only include the sample if it's not a fill value
                  if (this%num_samples(ind1) == 0) then
                     this%data(ind1) = fld_val
                     this%var_buffer(ind1) = 0.0_REAL64
                     this%num_samples(ind1) = this%num_samples(ind1) + 1
                  else
                     tmp = this%data(ind1)
                     this%num_samples(ind1) = this%num_samples(ind1) + 1
                     this%data(ind1) = this%data(ind1) + &
                             (fld_val - this%data(ind1)) / &
                             this%num_samples(ind1)
                     this%var_buffer(ind1) = this%var_buffer(ind1) + &
                             (fld_val - this%data(ind1)) * &
                             (fld_val - tmp)
                  end if
               end if
            else
               if (this%num_samples(ind1) == 0) then
                  this%data(ind1) = fld_val
                  this%var_buffer(ind1) = 0.0_REAL64
                  this%num_samples(ind1) = this%num_samples(ind1) + 1
               else
                  tmp = this%data(ind1)
                  this%num_samples(ind1) = this%num_samples(ind1) + 1
                  this%data(ind1) = this%data(ind1) + &
                          (fld_val - this%data(ind1)) / &
                          this%num_samples(ind1)
                  this%var_buffer(ind1) = this%var_buffer(ind1) + &
                          (fld_val - this%data(ind1)) * &
                          (fld_val - tmp)
               end if
            end if
         end do
      case default
         ! This shouldn't trigger, so raise an error:
         write(errmsg, *) 'un-recognized accumulation type: ', this%accum_type
         call hist_add_error(subname, errmsg, errors=logger)
      end select

   end subroutine buff_1d_accum

   !#######################################################################

   subroutine buff_1d_value(this, norm_val, flag_xyfill, fill_value, logger)
      use hist_msg_handler, only: hist_log_messages, ERROR, VERBOSE
      ! Dummy arguments
      class(hist_buff_1d_t),             intent(inout) :: this
      real(REAL64),                      intent(inout) :: norm_val(:)
      logical,                           intent(in)    :: flag_xyfill
      real(REAL64),                      intent(in)    :: fill_value
      type(hist_log_messages), optional, intent(inout) :: logger
      ! Local variable
      integer :: ind1
      integer :: nacc
      real(REAL64) :: tmpfill, variance

      if (this%accum_type /= hist_accum_var) then
         do ind1 = 1, this%field_shape(1)
            nacc = this%num_samples(ind1)
            if (nacc > 0) then
               norm_val(ind1) = this%data(ind1) / nacc
            else if (flag_xyfill) then
               norm_val(ind1) = fill_value
            end if
         end do
      else
         ! Standard deviation using Welford's algorithm
         ! DOI: 10.1080/00401706.1962.10490022
         tmpfill = merge(fill_value, 0.0_REAL64, flag_xyfill)
         do ind1 = 1, this%field_shape(1)
            if (this%num_samples(ind1) > 0) then
               variance = this%var_buffer(ind1) / this%num_samples(ind1)
               norm_val(ind1) = sqrt(variance)
            else
               norm_val(ind1) = tmpfill
            end if
         end do
      end if

   end subroutine buff_1d_value

   !#######################################################################

   subroutine buff_2d_clear(this, logger)
      use hist_msg_handler, only: hist_log_messages, hist_add_alloc_error

      ! Dummy arguments
      class(hist_buff_2d_t),       intent(inout) :: this
      type(hist_log_messages), optional, intent(inout) :: logger
      ! Local variables
      integer                     :: aerr
      character(len=256)          :: errmsg
      character(len=*), parameter :: subname = 'buff_2d_clear'

      if (.not. allocated(this%data)) then
         allocate(this%data(this%field_shape(1), this%field_shape(2)), stat=aerr, errmsg=errmsg)
         if (aerr /= 0) then
            call hist_add_alloc_error('data, error: '//errmsg, __FILE__, __LINE__ - 1, &
                 subname=subname, errors=logger)
         end if
      end if
      if (.not. allocated(this%var_buffer) .and. this%accum_type == hist_accum_var) then
         allocate(this%var_buffer(this%field_shape(1), this%field_shape(2)), stat=aerr, errmsg=errmsg)
         if (aerr /= 0) then
            call hist_add_alloc_error('var_buffer, error: '//errmsg, __FILE__, __LINE__ - 2,  &
                 subname=subname, errors=logger)
         end if
      end if
      if (.not. allocated(this%num_samples)) then
         allocate(this%num_samples(this%field_shape(1)), stat=aerr, errmsg=errmsg)
         if (aerr /= 0) then
            call hist_add_alloc_error('num_samples, error: '//errmsg, __FILE__, __LINE__ - 1, &
                 subname=subname, errors=logger)
         end if
      end if
      this%data = 0.0_REAL64
      this%num_samples = 0
      if (allocated(this%var_buffer)) then
         this%var_buffer = 0.0_REAL64
      end if

   end subroutine buff_2d_clear

   !#######################################################################

   subroutine init_buff_2d(this, field_in, volume_in, horiz_axis_in, &
        accum_type_in, shape_in, block_sizes_in, block_ind_in, logger)
      use hist_msg_handler, only: hist_log_messages

      class(hist_buff_2d_t),             intent(inout) :: this
      class(hist_hashable_t),  pointer,  intent(in)    :: field_in
      integer,                           intent(in)    :: volume_in
      integer,                           intent(in)    :: horiz_axis_in
      integer,                           intent(in)    :: accum_type_in
      integer,                           intent(in)    :: shape_in(:)
      integer,                 optional, intent(in)    :: block_sizes_in(:)
      integer,                 optional, intent(in)    :: block_ind_in
      type(hist_log_messages), optional, intent(inout) :: logger

      call init_buffer(this, field_in, volume_in, horiz_axis_in, accum_type_in, &
              shape_in, block_sizes_in, block_ind_in, logger=logger)
      call this%clear(logger=logger)
      this%buff_type = 'buff_2d'

   end subroutine init_buff_2d

   !#######################################################################

   subroutine buff_2d_accum(this, field, cols_or_block, flag_xyfill, &
         fill_value, cole, logger)
      use hist_msg_handler, only: hist_log_messages, hist_add_error
      ! Dummy arguments
      class(hist_buff_2d_t),       intent(inout) :: this
      real(REAL64),                      intent(in)    :: field(:,:)
      integer,                           intent(in)    :: cols_or_block
      logical,                           intent(in)    :: flag_xyfill
      real(REAL64),                      intent(in)    :: fill_value
      integer,                 optional, intent(in)    :: cole
      type(hist_log_messages), optional, intent(inout) :: logger
      ! Local variables
      integer            :: col_beg_use
      integer            :: col_end_use
      integer            :: ind1, ind2
      real(REAL64)       :: fld_val, tmp
      character(len=256) :: errmsg
      character(len=*), parameter :: subname = 'buff_2d_accum'

      if (this%has_blocks()) then
         ! For a blocked field, <cols_or_block> is a block index
         col_beg_use = this%block_begs(cols_or_block)
         col_end_use = this%block_ends(cols_or_block)
      else
         ! Non blocked, <cols_or_block> is the first column index
         col_beg_use = cols_or_block
         if (present(cole)) then
            col_end_use = cole
         else
            col_end_use = col_beg_use +                                       &
                 this%field_shape(this%horiz_axis_ind) - 1
         end if
      end if

      select case (this%accum_type)
      case (hist_accum_lst)
         do ind1 = col_beg_use, col_end_use
            ! Only set samples for this column if not a fill value
            if (flag_xyfill .and. field(ind1 - col_beg_use + 1, 1) == fill_value) then
               this%num_samples(ind1) = 0
            else
               this%num_samples(ind1) = 1
            end if
            do ind2 = 1, this%field_shape(2)
               fld_val = field(ind1 - col_beg_use + 1, ind2)
               if (flag_xyfill .and. fld_val == fill_value) then
                  this%data(ind1, ind2) = fill_value
               else
                  this%data(ind1, ind2) = fld_val
               end if
            end do
         end do
         if (flag_xyfill) then
            ! Check if assumption that columns are consistent wrt fill value was accurate
            call this%check_fill_value(field(col_beg_use:col_end_use, :), fill_value, logger)
         end if
      case (hist_accum_min)
         do ind1 = col_beg_use, col_end_use
            ! Has this column been sampled yet?
            if (this%num_samples(ind1) == 0) then
               ! First sample for this column: initialize all levels.
               ! Skip fill value columns (num_samples stays 0;
               ! buff_2d_value will write fill_value at output time)
               if (.not. (flag_xyfill .and. field(ind1 - col_beg_use + 1, 1) == fill_value)) then
                  do ind2 = 1, this%field_shape(2)
                     this%data(ind1, ind2) = field(ind1 - col_beg_use + 1, ind2)
                  end do
                  this%num_samples(ind1) = 1
               end if
            else
               ! Subsequent samples in column: update minimum per level
               do ind2 = 1, this%field_shape(2)
                  fld_val = field(ind1 - col_beg_use + 1, ind2)
                  if (flag_xyfill) then
                     if (fld_val < this%data(ind1, ind2) .and. &
                         fld_val /= fill_value) then
                        this%data(ind1, ind2) = fld_val
                     end if
                  else
                     if (fld_val < this%data(ind1, ind2)) then
                        this%data(ind1, ind2) = fld_val
                     end if
                  end if
               end do
            end if
         end do
         if (flag_xyfill) then
            call this%check_fill_value(field(col_beg_use:col_end_use, :), fill_value, logger)
         end if
      case (hist_accum_max)
         do ind1 = col_beg_use, col_end_use
            ! Has this column been sampled yet?
            if (this%num_samples(ind1) == 0) then
               ! First sample for this column: initialize all levels.
               ! Skip fill value columns (num_samples stays 0;
               ! buff_2d_value will write fill_value at output time)
               if (.not. (flag_xyfill .and. field(ind1 - col_beg_use + 1, 1) == fill_value)) then
                  do ind2 = 1, this%field_shape(2)
                     this%data(ind1, ind2) = field(ind1 - col_beg_use + 1, ind2)
                  end do
                  this%num_samples(ind1) = 1
               end if
            else
               ! Subsequent samples: update maximum per level
               do ind2 = 1, this%field_shape(2)
                  fld_val = field(ind1 - col_beg_use + 1, ind2)
                  if (flag_xyfill) then
                     if (fld_val > this%data(ind1, ind2) .and. fld_val /= fill_value) then
                        this%data(ind1, ind2) = fld_val
                     end if
                  else
                     if (fld_val > this%data(ind1, ind2)) then
                        this%data(ind1, ind2) = fld_val
                     end if
                  end if
               end do
            end if
         end do
         if (flag_xyfill) then
            call this%check_fill_value(field(col_beg_use:col_end_use, :), fill_value, logger)
         end if
      case (hist_accum_avg)
         do ind1 = col_beg_use, col_end_use
            if (.not. (flag_xyfill .and. field(ind1 - col_beg_use + 1, 1) == fill_value)) then
               this%num_samples(ind1) = this%num_samples(ind1) + 1
            end if
            do ind2 = 1, this%field_shape(2)
               fld_val = field(ind1 - col_beg_use + 1, ind2)
               ! Compute running sum
               if (flag_xyfill) then
                  ! Only include sample if it is not the fill value
                  if (fld_val /= fill_value) then
                     this%data(ind1, ind2) = this%data(ind1, ind2) + fld_val
                  end if
               else
                  this%data(ind1, ind2) = this%data(ind1, ind2) + fld_val
               end if
            end do
         end do
         if (flag_xyfill) then
            call this%check_fill_value(field(col_beg_use:col_end_use, :), fill_value, logger)
         end if
      case (hist_accum_var)
         ! Standard deviation using Welford's algorithm
         ! DOI: 10.1080/00401706.1962.10490022
         do ind1 = col_beg_use, col_end_use
            if (flag_xyfill .and. field(ind1 - col_beg_use + 1, 1) == fill_value) then
               this%num_samples(ind1) = this%num_samples(ind1)
            else
               this%num_samples(ind1) = this%num_samples(ind1) + 1
            end if
            do ind2 = 1, this%field_shape(2)
               fld_val = field(ind1 - col_beg_use + 1, ind2)
               if (flag_xyfill) then
                  if (fld_val /= fill_value) then
                     ! Only include the sample if it's not a fill value
                     if (this%num_samples(ind1) == 1) then
                        this%data(ind1, ind2) = fld_val
                        this%var_buffer(ind1, ind2) = 0.0_REAL64
                     else
                        tmp = this%data(ind1, ind2)
                        this%data(ind1, ind2) = this%data(ind1, ind2) + &
                                (fld_val - this%data(ind1, ind2)) / &
                                this%num_samples(ind1)
                        this%var_buffer(ind1, ind2) = this%var_buffer(ind1, ind2) + &
                                (fld_val - this%data(ind1, ind2)) * &
                                (fld_val - tmp)
                     end if
                  end if
               else
                  if (this%num_samples(ind1) == 1) then
                     this%data(ind1, ind2) = fld_val
                     this%var_buffer(ind1, ind2) = 0.0_REAL64
                  else
                     tmp = this%data(ind1, ind2)
                     this%data(ind1, ind2) = this%data(ind1, ind2) + &
                             (fld_val - this%data(ind1, ind2)) / &
                             this%num_samples(ind1)
                     this%var_buffer(ind1, ind2) = this%var_buffer(ind1, ind2) + &
                             (fld_val - this%data(ind1, ind2)) * &
                             (fld_val - tmp)
                  end if
               end if
            end do
         end do
         if (flag_xyfill) then
            call this%check_fill_value(field(col_beg_use:col_end_use, :), fill_value, logger)
         end if
      case default
         ! This shouldn't trigger, so raise an error:
         write(errmsg, *) 'un-recognized accumulation type: ', this%accum_type
         call hist_add_error(subname, errmsg, errors=logger)
      end select

   end subroutine buff_2d_accum

   !#######################################################################

   subroutine buff_2d_check_fill(this, field, fill_value, logger)
      use hist_msg_handler, only: hist_log_messages, hist_add_error, ERROR
      ! Ensure that columns have fill value applied consistently across levels
      class(hist_buff_2d_t),             intent(inout) :: this
      real(REAL64),                      intent(in)    :: field(:,:)
      real(REAL64),                      intent(in)    :: fill_value
      type(hist_log_messages), optional, intent(inout) :: logger
      ! Local variables
      character(len=512) :: errstr
      integer :: idx, jdx
      logical :: error_found

      error_found = .false.
      lev_loop: do jdx = 2, this%field_shape(2)
         col_loop: do idx = 1, this%field_shape(1)
            if (field(idx,1) == fill_value .and. field(idx,jdx) /= fill_value .or. &
                 field(idx,1) /= fill_value .and. field(idx,jdx) == fill_value) then
               write(errstr, '(a,i0)') 'ERROR: fill value applied inconsistently for column ', idx
               call hist_add_error('buff_2d_check_fill', errstr, errors=logger)
               error_found = .true.
               exit col_loop
            end if
         end do col_loop
         if (error_found) then
            exit lev_loop
         end if
      end do lev_loop

   end subroutine buff_2d_check_fill

   !#######################################################################

   subroutine buff_2d_value(this, norm_val, flag_xyfill, fill_value, logger)
      use hist_msg_handler, only: hist_log_messages, ERROR, VERBOSE
      ! Dummy arguments
      class(hist_buff_2d_t),             intent(inout) :: this
      real(REAL64),                      intent(inout) :: norm_val(:,:)
      logical,                           intent(in)    :: flag_xyfill
      real(REAL64),                      intent(in)    :: fill_value
      type(hist_log_messages), optional, intent(inout) :: logger
      ! Local variable
      integer :: ind1, ind2
      integer :: nacc
      real(REAL64) :: tmpfill, variance

      if (this%accum_type /= hist_accum_var) then
         do ind1 = 1, this%field_shape(1)
            do ind2 = 1, this%field_shape(2)
               nacc = this%num_samples(ind1)
               if (nacc > 0) then
                  norm_val(ind1, ind2) = this%data(ind1, ind2) / nacc
               else if (flag_xyfill) then
                  norm_val(ind1, ind2) = fill_value
               end if
            end do
         end do
      else
         ! Standard deviation using Welford's algorithm
         ! DOI: 10.1080/00401706.1962.10490022
         tmpfill = merge(fill_value, 0.0_REAL64, flag_xyfill)
         do ind1 = 1, this%field_shape(1)
            do ind2 = 1, this%field_shape(2)
               if (this%num_samples(ind1) > 0) then
                  variance = this%var_buffer(ind1, ind2) / this%num_samples(ind1)
                  norm_val(ind1, ind2) = sqrt(variance)
               else
                  norm_val(ind1, ind2) = tmpfill
               end if
            end do
         end do
      end if

   end subroutine buff_2d_value

   !#######################################################################

   function buffer_factory(buffer_type, logger) result(newbuf)
      ! Create a new buffer based on <buffer_type>.
      ! <buffer_type> has a format typekind_rank_accum
      ! Where:
      ! <typekind> is a lowercase string representation of a
      !    supported kind from the ISO_FORTRAN_ENV module.
      ! <rank> is the rank of the buffer (no leading zeros)
      ! <accum> is the accumulation type, one of:
      !    lst: Store the last value collected
      !    avg: Accumulate running average
      !    var: Accumulate standard deviation
      !    min: Accumulate smallest value
      !    max: Accumulate largest value

      use hist_msg_handler, only: hist_log_messages, ERROR, VERBOSE
      use hist_msg_handler, only: hist_add_error

      ! Arguments
      class(hist_buffer_t),    pointer                 :: newbuf
      character(len=*),                  intent(in)    :: buffer_type
      type(hist_log_messages), optional, intent(inout) :: logger

      ! Local variables
      character(len=*),   parameter :: subname = 'buffer_factory'
      character(len=512)            :: alloc_errmsg
      ! For buffer allocation
      integer                                    :: aerr
      type(hist_buff_1d_t), pointer :: real_1
      type(hist_buff_2d_t), pointer :: real_2

      ! Clear relevant pointers
      nullify(newbuf)
      nullify(real_1)
      nullify(real_2)

      ! Create new buffer
      select case (trim(buffer_type))
      case ('real_1')
         allocate(real_1, stat=aerr, errmsg=alloc_errmsg)
         if (aerr == 0) then
            newbuf => real_1
         else
            call hist_add_error(subname, &
               'Failed to allocate real_1 buffer, errmsg = ', &
               errstr2=trim(alloc_errmsg), errors=logger)
         end if
      case ('real_2')
         allocate(real_2, stat=aerr, errmsg=alloc_errmsg)
         if (aerr == 0) then
            newbuf => real_2
         else
            call hist_add_error(subname, &
               'Failed to allocate real_2 buffer, errmsg = ', &
               errstr2=trim(alloc_errmsg), errors=logger)
         end if
      case default
         call hist_add_error(subname,                                         &
              "Invalid or unsupported buffer type, '",                        &
              errstr2=trim(buffer_type), errstr3="'", errors=logger)
      end select
   end function buffer_factory

   !#######################################################################

end module hist_buffer
