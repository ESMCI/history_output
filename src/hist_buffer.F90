module hist_buffer
   use ISO_FORTRAN_ENV, only: REAL64, REAL32, INT32, INT64
   use hist_hashable,   only: hist_hashable_t

   implicit none
   private

   ! Public interfaces
   public :: buffer_factory

   ! Processing flag indices
   integer, parameter :: hist_proc_last    = 1 ! Save last sample
   integer, parameter :: hist_proc_average = 2 ! Average samples
   integer, parameter :: hist_proc_stddev  = 3 ! Standard deviation of samples
   integer, parameter :: hist_proc_min     = 4 ! Minimum of samples
   integer, parameter :: hist_proc_max     = 5 ! Maximum of samples
   ! Time sampling flag indices
   !!XXgoldyXX: Todo: decide on sampling types

   type, abstract, public :: hist_buffer_t
      ! hist_buffer_t is an abstract base class for hist_outfld buffers
      class(hist_hashable_t), pointer              :: field_info => NULL()
      integer,                             private :: vol = -1 ! For host output
      integer,                             private :: horiz_axis_ind = 0
      character(len=:),       allocatable, private :: accum_str
      class(hist_buffer_t),   pointer              :: next
   contains
      procedure                                 :: field  => get_field_info
      procedure                                 :: volume => get_volume
      procedure                                 :: horiz_axis_index
      procedure                                 :: init_buffer
      procedure                                 :: accum_string
      procedure(hist_buff_init),       deferred :: initialize
      procedure(hist_buff_sub_noargs), deferred :: clear
   end type hist_buffer_t

   type, public, extends(hist_buffer_t) :: hist_buffer_1dreal32_inst_t
      integer               :: num_samples = 0
      real(REAL32), pointer :: data(:) => NULL()
   CONTAINS
      procedure :: clear => buff_1dreal32_inst_clear
      procedure :: accumulate => buff_1dreal32_inst_accum
      procedure :: norm_value => buff_1dreal32_inst_value
      procedure :: initialize => init_buff_1dreal32
   end type hist_buffer_1dreal32_inst_t

   type, public, extends(hist_buffer_t) :: hist_buffer_1dreal64_inst_t
      integer               :: num_samples = 0
      real(REAL64), pointer :: data(:) => NULL()
   CONTAINS
      procedure :: clear => buff_1dreal64_inst_clear
      procedure :: accumulate => buff_1dreal64_inst_accum
      procedure :: norm_value => buff_1dreal64_inst_value
      procedure :: initialize => init_buff_1dreal64
   end type hist_buffer_1dreal64_inst_t

   ! Abstract interfaces for hist_buffer_t class
   abstract interface
      subroutine hist_buff_sub_noargs(this)
         import :: hist_buffer_t
         class(hist_buffer_t), intent(inout) :: this
      end subroutine hist_buff_sub_noargs
   end interface

   abstract interface
      subroutine hist_buff_init(this, field_in, volume_in, horiz_axis_in,     &
           shape_in, block_sizes_in, block_ind_in)
         import :: hist_buffer_t
         import :: hist_hashable_t
         class(hist_buffer_t), intent(inout) :: this
         class(hist_hashable_t), pointer     :: field_in
         integer,              intent(in)    :: volume_in
         integer,              intent(in)    :: horiz_axis_in
         integer,              intent(in)    :: shape_in(:)
         integer,              intent(in)    :: block_sizes_in(:)
         integer,              intent(in)    :: block_ind_in
      end subroutine hist_buff_init
   end interface

CONTAINS

   function get_field_info(this)
      class(hist_buffer_t), intent(in) :: this
      class(hist_hashable_t), pointer  :: get_field_info

      get_field_info => this%field_info
   end function get_field_info

   !#######################################################################

   integer function get_volume(this)
      class(hist_buffer_t), intent(in) :: this

      get_volume = this%vol
   end function get_volume

   !#######################################################################

   integer function horiz_axis_index(this)
      class(hist_buffer_t), intent(in) :: this

      horiz_axis_index = this%horiz_axis_ind
   end function horiz_axis_index

   !#######################################################################

   subroutine init_buffer(this, field_in, volume_in, horiz_axis_in,           &
        block_sizes_in, block_ind_in)
      class(hist_buffer_t), intent(inout) :: this
      class(hist_hashable_t), pointer     :: field_in
      integer,              intent(in)    :: volume_in
      integer,              intent(in)    :: horiz_axis_in
      integer,              intent(in)    :: block_sizes_in(:)
      integer,              intent(in)    :: block_ind_in

      this%field_info => field_in
      this%vol = volume_in
      this%horiz_axis_ind = horiz_axis_in
   end subroutine init_buffer

   !#######################################################################

   function accum_string(this) result(ac_str)
      class(hist_buffer_t), intent(in) :: this
      character(len=:), allocatable    :: ac_str

      ac_str = this%accum_str
   end function accum_string

   !#######################################################################

   subroutine buff_1dreal32_inst_clear(this)
      class(hist_buffer_1dreal32_inst_t), intent(inout) :: this

      this%num_samples = 0
   end subroutine buff_1dreal32_inst_clear

   !#######################################################################

   subroutine init_buff_1dreal32(this, field_in, volume_in, horiz_axis_in, &
        shape_in, block_sizes_in, block_ind_in)
      class(hist_buffer_1dreal32_inst_t), intent(inout) :: this
      class(hist_hashable_t), pointer                   :: field_in
      integer,                            intent(in)    :: volume_in
      integer,                            intent(in)    :: horiz_axis_in
      integer,                            intent(in)    :: shape_in(:)
      integer,                            intent(in)    :: block_sizes_in(:)
      integer,                            intent(in)    :: block_ind_in

      call init_buffer(this, field_in, volume_in, horiz_axis_in,              &
           block_sizes_in, block_ind_in)
      this%accum_str = 'last sampled value'
      allocate(this%data(shape_in(1)))

   end subroutine init_buff_1dreal32

   !#######################################################################

   subroutine buff_1dreal32_inst_accum(this, field, errmsg)
      class(hist_buffer_1dreal32_inst_t), intent(inout) :: this
      real(REAL32),                       intent(in)    :: field(:)
      character(len=*), optional,         intent(out)   :: errmsg

      if (present(errmsg)) then
         errmsg = 'Not implemented'
      end if
      this%num_samples = 1

   end subroutine buff_1dreal32_inst_accum

   !#######################################################################

   subroutine buff_1dreal32_inst_value(this, norm_val, errmsg)
      class(hist_buffer_1dreal32_inst_t), intent(inout) :: this
      real(REAL32),                       intent(inout) :: norm_val(:)
      character(len=*), optional,         intent(out)   :: errmsg

      if (present(errmsg)) then
         errmsg = 'Not implemented'
      end if

   end subroutine buff_1dreal32_inst_value

   !#######################################################################

   subroutine buff_1dreal64_inst_clear(this)
      class(hist_buffer_1dreal64_inst_t), intent(inout) :: this

      this%num_samples = 0

   end subroutine buff_1dreal64_inst_clear

   !#######################################################################

   subroutine init_buff_1dreal64(this, field_in, volume_in, horiz_axis_in,    &
        shape_in, block_sizes_in, block_ind_in)
      class(hist_buffer_1dreal64_inst_t), intent(inout) :: this
      class(hist_hashable_t), pointer                   :: field_in
      integer,                            intent(in)    :: volume_in
      integer,                            intent(in)    :: horiz_axis_in
      integer,                            intent(in)    :: shape_in(:)
      integer,                            intent(in)    :: block_sizes_in(:)
      integer,                            intent(in)    :: block_ind_in

      call init_buffer(this, field_in, volume_in, horiz_axis_in,              &
           block_sizes_in, block_ind_in)
      this%accum_str = 'last sampled value'
      allocate(this%data(shape_in(1)))

   end subroutine init_buff_1dreal64

   !#######################################################################

   subroutine buff_1dreal64_inst_accum(this, field, errmsg)
      class(hist_buffer_1dreal64_inst_t), intent(inout) :: this
      real(REAL64),                       intent(in)    :: field(:)
      character(len=*), optional,         intent(out)   :: errmsg

      if (present(errmsg)) then
         errmsg = 'Not implemented'
      end if

   end subroutine buff_1dreal64_inst_accum

   !#######################################################################

   subroutine buff_1dreal64_inst_value(this, norm_val, errmsg)
      class(hist_buffer_1dreal64_inst_t), intent(inout) :: this
      real(REAL64),                       intent(inout) :: norm_val(:)
      character(len=*), optional,         intent(out)   :: errmsg

      if (present(errmsg)) then
         errmsg = 'Not implemented'
      end if

   end subroutine buff_1dreal64_inst_value

   !#######################################################################

   function buffer_factory(buffer_type) result(newbuf)
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
      ! Arguments
      class(hist_buffer_t), pointer :: newbuf
      character(len=*), intent(in)  :: buffer_type
      ! Local variables
      ! For buffer allocation
      type(hist_buffer_1dreal32_inst_t), pointer :: real32_1_in => NULL()
      type(hist_buffer_1dreal64_inst_t), pointer :: real64_1_in => NULL()


      ! Create new buffer
      select case (trim(buffer_type))
      case ('real32_1_inst')
         allocate(real32_1_in)
         newbuf => real32_1_in
      case ('real64_1_inst')
         allocate(real64_1_in)
         newbuf => real64_1_in
      end select
   end function buffer_factory

end module hist_buffer