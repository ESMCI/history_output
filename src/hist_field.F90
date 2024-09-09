module hist_field
   ! Module containing DDTs for history fields and associated routines

   use hist_hashable, only: hist_hashable_t
   use hist_buffer,   only: hist_buffer_t

   implicit none
   private

   public :: hist_field_initialize
   public :: hist_get_field ! Get field from a buffer

   type, public, extends(hist_hashable_t) :: hist_field_info_t
      ! Field metadata (note: all strings should be lowercase)
      character(len=:), allocatable, private :: diag_file_name
      character(len=:), allocatable, private :: field_standard_name
      character(len=:), allocatable, private :: field_long_name
      character(len=:), allocatable, private :: field_units
      character(len=:), allocatable, private :: field_type
      character(len=:), allocatable, private :: field_kind
      character(len=:), allocatable, private :: field_accumulate_type
      character(len=:), allocatable, private :: field_sampling_sequence
      character(len=:), allocatable, private :: field_mixing_ratio
      character(len=:), allocatable, private :: field_cell_methods
      integer,          allocatable, private :: field_shape(:)
      integer,                       private :: field_decomp
      integer,                       private :: field_num_levels
      integer,          allocatable, private :: field_dimensions(:)
      logical,                       private :: field_flag_xyfill
      integer,          allocatable, private :: field_beg_dims(:)
      integer,          allocatable, private :: field_end_dims(:)
      
! type kind rank?
      ! dimensions?
      type(hist_field_info_t), pointer :: next => NULL()
      class(hist_buffer_t),    pointer :: buffers => NULL()
   contains
      procedure :: key             => hist_field_info_get_key
      procedure :: diag_name       => get_diag_name
      procedure :: standard_name   => get_standard_name
      procedure :: long_name       => get_long_name
      procedure :: units           => get_units
      procedure :: type            => get_type
      procedure :: kind            => get_kind
      procedure :: shape           => get_shape
      procedure :: num_levels      => get_num_levels
      procedure :: accumulate_type => get_accumulate_type
      procedure :: decomp          => get_decomp
      procedure :: dimensions      => get_dimensions
      procedure :: beg_dims        => get_beg_dims
      procedure :: end_dims        => get_end_dims
      procedure :: sampling_sequence => get_sampling_sequence
      procedure :: flag_xyfill       => get_flag_xyfill
      procedure :: mixing_ratio      => get_mixing_ratio
      procedure :: cell_methods      => get_cell_methods
      procedure :: set_dimension_bounds
      procedure :: clear_buffers
      final     :: finalize_field
   end type hist_field_info_t

CONTAINS

   !#######################################################################

   function hist_field_info_get_key(hashable)
      ! Return the hashable field info class key (diag_file_name)
      class(hist_field_info_t), intent(in) :: hashable
      character(len=:), allocatable        :: hist_field_info_get_key

      hist_field_info_get_key = hashable%diag_file_name
   end function hist_field_info_get_key

   !#######################################################################

   subroutine hist_field_initialize(field, diag_name_in, std_name_in,         &
        long_name_in, units_in, type_in, decomp_in, mdim_indices, acc_type, num_levels,  &
        field_shape, sampling_seq, flag_xyfill, mixing_ratio, dim_bounds, mdim_sizes, &
        beg_dims, end_dims, cell_methods, errmsg)

      type(hist_field_info_t), pointer               :: field
      character(len=*),                  intent(in)  :: diag_name_in
      character(len=*),                  intent(in)  :: std_name_in
      character(len=*),                  intent(in)  :: long_name_in
      character(len=*),                  intent(in)  :: units_in
      character(len=*),                  intent(in)  :: type_in
      integer,                           intent(in)  :: decomp_in
      integer,                           intent(in)  :: mdim_indices(:)
      character(len=*),                  intent(in)  :: acc_type
      integer,                           intent(in)  :: num_levels
      integer,                           intent(in)  :: field_shape(:)
      character(len=*),        optional, intent(in)  :: sampling_seq
      logical,                 optional, intent(in)  :: flag_xyfill
      character(len=*),        optional, intent(in)  :: mixing_ratio
      integer,                 optional, intent(in)  :: dim_bounds(:,:)
      integer,                 optional, intent(in)  :: mdim_sizes(:)
      integer,                 optional, intent(in)  :: beg_dims(:)
      integer,                 optional, intent(in)  :: end_dims(:)
      character(len=*),        optional, intent(in)  :: cell_methods
      character(len=*),        optional, intent(out) :: errmsg

      if (present(errmsg)) then
         errmsg = ''
      end if
      field%diag_file_name = diag_name_in
      field%field_standard_name = std_name_in
      field%field_long_name = long_name_in
      field%field_units = units_in
      field%field_type = type_in
      field%field_decomp = decomp_in
      field%field_accumulate_type = acc_type
      field%field_num_levels = num_levels
      allocate(field%field_dimensions(size(mdim_indices, 1)))
      field%field_dimensions = mdim_indices
      allocate(field%field_shape(size(field_shape, 1)))
      field%field_shape = field_shape
      if (present(sampling_seq)) then
         field%field_sampling_sequence = sampling_seq
      else
         field%field_sampling_sequence = ''
      end if
      if (present(flag_xyfill)) then
         field%field_flag_xyfill = flag_xyfill
      else
         field%field_flag_xyfill = .false.
      end if
      if (present(mixing_ratio)) then
         field%field_mixing_ratio = mixing_ratio
      else
         field%field_mixing_ratio = ''
      end if
      if (present(cell_methods)) then
         field%field_cell_methods = cell_methods
      else
         field%field_cell_methods = ''
      end if
      if (present(dim_bounds)) then
         if (.not. present(mdim_sizes)) then
            if (present(errmsg)) then
               errmsg = 'hist_field_initialize: mdim_sizes must be provided for dimension bound setting'
               return
            end if
         end if
         call field%set_dimension_bounds(dim_bounds, mdim_sizes)
      end if
      if (present(beg_dims) .or. present(end_dims)) then
         if (.not. (present(beg_dims) .and. present(end_dims))) then
            errmsg = 'hist_field_initialize: both beg_dims and end_dims must be provided if one is'
            return
         end if
         if (.not. present(dim_bounds)) then
            allocate(field%field_beg_dims(size(beg_dims)))
            field%field_beg_dims = beg_dims
            allocate(field%field_end_dims(size(end_dims)))
            field%field_end_dims = end_dims
         end if
      end if
   end subroutine hist_field_initialize

   !#######################################################################

   subroutine hist_get_field(buffer, field)
      ! Retrieve
      class(hist_buffer_t), intent(in) :: buffer
      type(hist_field_info_t), pointer :: field

      nullify(field)
      select type(finfo => buffer%field_info)
      type is (hist_field_info_t)
         field => finfo
      end select
   end subroutine hist_get_field

   !#######################################################################

   function get_diag_name(this) result(info)
      class(hist_field_info_t), intent(in) :: this
      character(len=:), allocatable        :: info

      info = this%diag_file_name
   end function get_diag_name

   !#######################################################################

   function get_standard_name(this) result(info)
      class(hist_field_info_t), intent(in) :: this
      character(len=:), allocatable        :: info

      info = this%field_standard_name
   end function get_standard_name

   !#######################################################################

   function get_long_name(this) result(info)
      class(hist_field_info_t), intent(in) :: this
      character(len=:), allocatable        :: info

      info = this%field_long_name
   end function get_long_name

   !#######################################################################

   function get_units(this) result(info)
      class(hist_field_info_t), intent(in) :: this
      character(len=:), allocatable        :: info

      info = this%field_units
   end function get_units

   !#######################################################################

   function get_type(this) result(info)
      class(hist_field_info_t), intent(in) :: this
      character(len=:), allocatable        :: info

      info = this%field_type
   end function get_type

   !#######################################################################

   function get_kind(this) result(info)
      class(hist_field_info_t), intent(in) :: this
      character(len=:), allocatable        :: info

      info = this%field_kind
   end function get_kind

   !#######################################################################

   function get_accumulate_type(this) result(info)
      class(hist_field_info_t), intent(in) :: this
      character(len=:), allocatable        :: info

      info = this%field_accumulate_type
   end function get_accumulate_type

   !#######################################################################

   function get_decomp(this) result(decomp)
      class(hist_field_info_t), intent(in) :: this
      integer                              :: decomp

      decomp = this%field_decomp
   end function get_decomp

   !#######################################################################

   function get_num_levels(this) result(num_levels)
      class(hist_field_info_t), intent(in) :: this
      integer                              :: num_levels

      num_levels = this%field_num_levels
   end function get_num_levels

   !#######################################################################

   function get_dimensions(this) result(dimensions)
      class(hist_field_info_t), intent(in) :: this
      integer,   allocatable               :: dimensions(:)
      allocate(dimensions(size(this%field_dimensions)))
      dimensions = this%field_dimensions
   end function get_dimensions

   !#######################################################################

   function get_beg_dims(this) result(beg_dim)
      class(hist_field_info_t), intent(in) :: this
      integer,   allocatable               :: beg_dim(:)
      if (allocated(this%field_beg_dims)) then
         allocate(beg_dim(size(this%field_beg_dims)))
         beg_dim = this%field_beg_dims
      end if
   end function get_beg_dims

   !#######################################################################

   function get_end_dims(this) result(end_dim)
      class(hist_field_info_t), intent(in) :: this
      integer,   allocatable               :: end_dim(:)
      if (allocated(this%field_end_dims)) then
         allocate(end_dim(size(this%field_end_dims)))
         end_dim = this%field_end_dims
      end if
   end function get_end_dims

   !#######################################################################

   function get_shape(this) result(field_shape)
      class(hist_field_info_t), intent(in) :: this
      integer,   allocatable               :: field_shape(:)
      allocate(field_shape(size(this%field_shape)))
      field_shape = this%field_shape
   end function get_shape

   !#######################################################################

   function get_sampling_sequence(this) result(info)
      class(hist_field_info_t), intent(in) :: this
      character(len=:), allocatable        :: info

      info = this%field_sampling_sequence
   end function get_sampling_sequence

   !#######################################################################

   function get_flag_xyfill(this) result(flag_xyfill)
      class(hist_field_info_t), intent(in) :: this
      logical                              :: flag_xyfill

      flag_xyfill = this%field_flag_xyfill

   end function get_flag_xyfill

   !#######################################################################

   function get_mixing_ratio(this) result(info)
      class(hist_field_info_t), intent(in) :: this
      character(len=:), allocatable        :: info

      info = this%field_mixing_ratio
   end function get_mixing_ratio

   !#######################################################################

   function get_cell_methods(this) result(info)
      class(hist_field_info_t), intent(in) :: this
      character(len=:), allocatable        :: info

      info = this%field_cell_methods
   end function get_cell_methods

   !#######################################################################

   subroutine set_dimension_bounds(this, dimbounds, mdim_sizes)
      class(hist_field_info_t), intent(inout) :: this
      integer, intent(in) :: dimbounds(:,:)
      integer, intent(in) :: mdim_sizes(:)

      integer ::idx

      allocate(this%field_beg_dims(3))
      allocate(this%field_end_dims(3))

      this%field_beg_dims(1) = dimbounds(1,1)
      this%field_end_dims(1) = dimbounds(1,2)
      this%field_beg_dims(2) = 1

      if (allocated(this%field_dimensions)) then
         if (size(this%field_dimensions) > 0) then
            this%field_end_dims(2) = 1
            do idx = 1, size(this%field_dimensions)
               this%field_end_dims(2) = this%field_end_dims(2) * mdim_sizes(idx)
            end do
         else
            this%field_end_dims(2) = this%field_num_levels
         end if
      else
         this%field_end_dims(2) = this%field_num_levels
      end if

      this%field_beg_dims(3) = dimbounds(2,1)
      this%field_end_dims(3) = dimbounds(2,2)

   end subroutine set_dimension_bounds

   !#######################################################################

   subroutine clear_buffers(this, logger)
      use hist_msg_handler,    only: hist_log_messages
      ! Dummy Argument
      class(hist_field_info_t) :: this
      type(hist_log_messages), optional :: logger
      ! Local Variables
      class(hist_buffer_t), pointer :: next_buf

      next_buf => this%buffers
      do
         if (associated(next_buf)) then
            call next_buf%clear(logger=logger)
            next_buf => next_buf%next
         else
            exit
         end if
      end do

   end subroutine

   !#######################################################################

   subroutine finalize_field(this)
      ! Dummy Argument
      type(hist_field_info_t) :: this
      ! Local Variables
      class(hist_buffer_t), pointer :: next_buf

      if (allocated(this%diag_file_name)) then
         deallocate(this%diag_file_name)
      end if
      if (allocated(this%field_standard_name)) then
         deallocate(this%field_standard_name)
      end if
      if (allocated(this%field_long_name)) then
         deallocate(this%field_long_name)
      end if
      if (allocated(this%field_units)) then
         deallocate(this%field_units)
      end if
      if (allocated(this%field_type)) then
         deallocate(this%field_type)
      end if
      if (allocated(this%field_accumulate_type)) then
         deallocate(this%field_accumulate_type)
      end if
      if (allocated(this%field_dimensions)) then
         deallocate(this%field_dimensions)
      end if
      if (allocated(this%field_sampling_sequence)) then
         deallocate(this%field_sampling_sequence)
      end if
      if (allocated(this%field_shape)) then
         deallocate(this%field_shape)
      end if
      ! We are not in charge of the field chain so just nullify
      nullify(this%next)
      ! We are in charge of the buffers so get rid of them.
      do
         next_buf => this%buffers
         if (associated(next_buf)) then
            this%buffers => next_buf%next
            deallocate(next_buf)
         else
            exit
         end if
      end do
   end subroutine finalize_field

end module hist_field
