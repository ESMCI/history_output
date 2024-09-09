!!XXgoldyXX: To do, statistics output
module hist_hash_table

   use hist_hashable, only: hist_hashable_t

   implicit none
   private

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   !  Hashing.
   !
   !  Accelerate processing (e.g., of outfld) by using a hash function of
   !    the field name
   !
   !  Note: the hashing logic will fail if any of the following are true:
   !
   !         1) The lower bound on the dimension of 'masterlist' is less than 1.
   !
   !         2) 'outfld' is called with field names that are not defined on
   !            masterlist.  This applies to both initial/branch and restart
   !            runs.
   !
   !         3) An inconsistency between a field's file active flag
   !            'masterlist(ff)%actflag(t)' and active fields read from
   !            restart files.
   !
   !         4) Invoking function 'gen_hash_key' before the primary and
   !            secondary hash tables have been created
   !            (routine bld_outfld_hash_tbls).
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !
   !  User definable constants for hash and overflow tables.
   !  Define size of primary hash table (specified as 2**size).
   !
!   integer, parameter :: tbl_hash_pri_sz_lg2 = 16
   !
   !  Define size of overflow hash table % of primary hash table.
   !
!   integer, parameter :: tbl_hash_oflow_percent = 20
   !
   !  Do *not* modify the parameters below.
   !
!   integer, parameter :: tbl_hash_pri_sz = 2**tbl_hash_pri_sz_lg2
!   integer, parameter :: tbl_hash_oflow_sz = tbl_hash_pri_sz * (tbl_hash_oflow_percent/100.0_r8)
   !
   !  The primary and overflow tables are organized to mimimize space (read:
   !  try to maximimze cache line usage).
   !
   !  gen_hash_key(fieldname) will return an index on the interval
   !  [0 ... tbl_hash_pri_sz-1].
   !
   !
   !  Primary:
   !  gen_hash_key(fieldname)-------+     +----------+
   !                                |     |   -ii    | 1 ------>tbl_hash_oflow(ii)
   !                                |     +----------+
   !                                +-->  |    ff    | 2 ------>masterlist(ff)
   !                                      +----------+
   !                                      |          | ...
   !                                      +----------+
   !                                      |          | tbl_hash_pri_sz
   !                                      +----------+
   !
   !  Overflow (if tbl_hash_pri() < 0):
   !  tbl_hash_pri(gen_hash_key(fieldname))
   !                         |
   !                         |            +----------+
   !                         |            |     1    | 1  (one entry on O.F. chain)
   !                         |            +----------+
   !                         |            |    ff_m  | 2
   !                         |            +----------+
   !                         +--------->  |     3    | 3  (three entries on chain)
   !                                      +----------+
   !                                      |    ff_x  | 4
   !                                      +----------+
   !                                      |    ff_y  | 5
   !                                      +----------+
   !                                      |    ff_z  | 6
   !                                      +----------+
   !                                      |          | ...
   !                                      +----------+
   !                                      |          | tbl_hash_oflow_sz
   !                                      +----------+
   !
   !
   !
   !  Constants used in hashing function gen_hash_key.
   !  Note: if the constants in table 'tbl_gen_hash_key' below are modified,
   !        changes are required to routine 'gen_hash_key' because of specific
   !        logic in the routine that optimizes character strings of length 8.
   !

   integer, parameter :: gen_hash_key_offset = 21467 ! z'000053db'

   integer, parameter :: tbl_max_idx = 15
   integer, parameter, dimension(0:tbl_max_idx) :: tbl_gen_hash_key =         &
        (/ 61, 59, 53, 47, 43, 41, 37, 31, 29, 23, 17, 13, 11, 7, 3, 1 /)

   integer, parameter :: table_factor_size = 8     ! Table size / # entries
   integer, parameter :: table_overflow_factor = 4 ! # entries / Overflow size

   type :: table_entry_t
      ! Any table entry contains a key and a value
      class(hist_hashable_t), pointer             :: entry_value => NULL()
      type(table_entry_t),    pointer             :: next => NULL()
   contains
      procedure :: finalize_table_entry
   end type table_entry_t

   type, public :: hist_hash_table_t
      ! hist_hash_table_t contains all information to build and use a hash table
      ! It also keeps track of statistics such as collision frequency and size
      integer, private                          :: table_size = -1
      integer, private                          :: overflow_size = -1
      integer, private                          :: key_offset = gen_hash_key_offset
      type(table_entry_t), private, allocatable :: primary_table(:)
      ! Statistics
      integer, private                          :: num_keys = 0
      integer, private                          :: num_key_collisions = 0
      integer, private                          :: max_collision = 0
   contains
      procedure :: initialize   => hash_table_initialize_table
      procedure :: key_hash     => hash_table_key_hash
      procedure :: add_hash_key => hash_table_add_hash_key
      procedure :: table_value  => hash_table_table_value
      procedure :: deallocate_table => hash_table_deallocate_table
   end type hist_hash_table_t

   !! Private interfaces
   private :: have_error      ! Has a called routine detected an error?
   private :: clear_optstring ! Clear a string, if present

CONTAINS

   !#######################################################################

   logical function have_error(errmsg)
      ! Return .true. iff <errmsg> is present and contains text
      character(len=*), optional, intent(in) :: errmsg

      have_error = present(errmsg)
      if (have_error) then
         have_error = len_trim(errmsg) > 0
      end if
   end function have_error

   !#######################################################################

   subroutine clear_optstring(str)
      ! clear <str> if it is present
      character(len=*), optional, intent(inout) :: str

      if (present(str)) then
         str = ''
      end if
   end subroutine clear_optstring

   !#######################################################################

   subroutine finalize_table_entry(this)
      class(table_entry_t), intent(inout)  :: this

      nullify(this%entry_value) ! We do not own the memory
      if (associated(this%next)) then
         deallocate(this%next) ! This should invoke finalize recursively
         nullify(this%next)
      end if

   end subroutine finalize_table_entry

   !#######################################################################

   subroutine hash_table_initialize_table(this, tbl_size, key_off)
      ! Initialize this table.
      ! Dummy arguments
      class(hist_hash_table_t)      :: this
      integer,           intent(in) :: tbl_size   ! new table size
      integer, optional, intent(in) :: key_off    ! key offset

      ! Clear this table so it can be initialized
      if (allocated(this%primary_table)) then
         call this%deallocate_table()
      end if
      ! Avoid too-large tables
      this%table_size = ishft(1, MIN(tbl_size, 14))

      allocate(this%primary_table(this%table_size))
      if (present(key_off)) then
         this%key_offset = key_off
      end if
   end subroutine hash_table_initialize_table

   !#######################################################################

   subroutine hash_table_deallocate_table(this)
      ! Deallocate the hash table (if present)
      ! Dummy argument
      class(hist_hash_table_t)     :: this

      ! Local variable
      integer :: entry_index

      if (allocated(this%primary_table)) then
         do entry_index = 1, size(this%primary_table)
            call this%primary_table(entry_index)%finalize_table_entry()
         end do
         deallocate(this%primary_table)
      end if

   end subroutine hash_table_deallocate_table

   !#######################################################################

   integer function hash_table_key_hash(this, string, errmsg) result(hash_key)
      !
      !-----------------------------------------------------------------------
      !
      ! Purpose: Generate a hash key on the interval [0 .. tbl_hash_pri_sz-1]
      !          given a character string.
      !
      ! Algorithm is a variant of perl's internal hashing function.
      !
      !-----------------------------------------------------------------------
      !
      !
      !  Arguments:
      !
      class(hist_hash_table_t)                :: this
      character(len=*),           intent(in)  :: string
      character(len=*), optional, intent(out) :: errmsg
      character(len=*), parameter             :: subname = 'HASH_TABLE_KEY_HASH'
      !
      !  Local.
      !
      integer :: hash
      integer :: index
      integer :: ind_fact
      integer :: hash_fact

      hash = this%key_offset
      ind_fact = 0
      do index = 1, len_trim(string)
         ind_fact = ind_fact + 1
         if (ind_fact > tbl_max_idx) then
            ind_fact = 1
         end if
         hash_fact = tbl_gen_hash_key(ind_fact)
         hash = ieor(hash, (ichar(string(index:index)) * hash_fact))
      end do

      hash_key = iand(hash, this%table_size - 1) + 1
      if ((hash_key < 1) .or. (hash_key > this%table_size)) then
         if (present(errmsg)) then
            write(errmsg, '(2a,2(i0,a))') subname, ' ERROR: Key Hash, ',      &
                 hash_key, ' out of bounds, [1, ', this%table_size, ']'
         else
            write(6, '(2a,2(i0,a))') subname, ' ERROR: Key Hash, ',           &
                 hash_key, ' out of bounds, [1, ', this%table_size, ']'
            STOP 1
         end if
      end if

   end function hash_table_key_hash

   !#######################################################################

   function hash_table_table_value(this, key, errmsg) result(tbl_val)
      !
      !-----------------------------------------------------------------------
      !
      ! Purpose: Return the the key value of <key>
      !
      !          If the object is not found, return NULL
      !
      !-----------------------------------------------------------------------
      !
      !  Arguments:
      !
      class(hist_hash_table_t)                      :: this
      character(len=*),                 intent(in)  :: key
      character(len=*),       optional, intent(out) :: errmsg
      class(hist_hashable_t), pointer               :: tbl_val
      !
      !  Local.
      !
      integer                        :: hash_key
      type(table_entry_t), pointer   :: next_ptr
      character(len=*),    parameter :: subname = 'HASH_TABLE_TABLE_INDEX'

      call clear_optstring(errmsg)
      nullify(tbl_val)
      hash_key = this%key_hash(key, errmsg=errmsg)
      ASSOCIATE(tbl_entry => this%primary_table(hash_key))
         if (have_error(errmsg)) then
            errmsg = trim(errmsg)//', called from '//subname
         else if (associated(tbl_entry%entry_value)) then
            if (tbl_entry%entry_value%key() == trim(key)) then
               tbl_val => tbl_entry%entry_value
            else
               next_ptr => tbl_entry%next
               do
                  if (associated(next_ptr)) then
                     if (associated(next_ptr%entry_value)) then
                        if (next_ptr%entry_value%key() == trim(key)) then
                           tbl_val => next_ptr%entry_value
                           exit
                        end if
                     end if
                     next_ptr => next_ptr%next
                  else
                     exit
                  end if
               end do
            end if
         end if
      END ASSOCIATE

      if ((.not. associated(tbl_val)) .and. present(errmsg)) then
         if (.not. have_error(errmsg)) then ! Still need to test for empty
            write(errmsg, *) subname, ": No entry for '", trim(key), "'"
         end if
      end if

   end function hash_table_table_value

   !#######################################################################

   subroutine hash_table_add_hash_key(this, newval, errmsg)
      !
      !-----------------------------------------------------------------------
      !
      ! Purpose: Add <newval> to this hash table using its key
      !          Its key must not be an empty string
      !          It is an error to try to add a key more than once
      !
      !
      !-----------------------------------------------------------------------

      !  Dummy arguments:
      class(hist_hash_table_t)                      :: this
      class(hist_hashable_t), target                :: newval
      character(len=*),       optional, intent(out) :: errmsg
      ! Local variables
      integer                          :: hash_ind
      integer                          :: ovflw_len
      character(len=:),    allocatable :: newkey
      type(table_entry_t), pointer     :: next_ptr
      type(table_entry_t), pointer     :: new_entry
      character(len=*),    parameter   :: subname = 'HASH_TABLE_ADD_HASH_KEY'

      call clear_optstring(errmsg)
      newkey = newval%key()
      hash_ind = this%key_hash(newkey, errmsg=errmsg)
      ! Check for this entry
      if (have_error(errmsg)) then
         errmsg = trim(errmsg)//', called from '//subname
      else if (associated(this%table_value(newkey))) then
         if (present(errmsg)) then
            write(errmsg, *) subname, " ERROR: key, '", newkey,               &
                 "' already in table"
         end if
      else
         ASSOCIATE(tbl_entry => this%primary_table(hash_ind))
            if (associated(tbl_entry%entry_value)) then
               ! We have a collision, make a new entry
               allocate(new_entry)
               new_entry%entry_value => newval
               ! Now, find a spot
               if (associated(tbl_entry%next)) then
                  ovflw_len = 1
                  next_ptr => tbl_entry%next
                  do
                     if (associated(next_ptr%next)) then
                        ovflw_len = ovflw_len + 1
                        next_ptr => next_ptr%next
                     else
                        exit
                     end if
                  end do
                  ovflw_len = ovflw_len + 1
                  next_ptr%next => new_entry
               else
                  this%num_key_collisions = this%num_key_collisions + 1
                  tbl_entry%next => new_entry
                  ovflw_len = 1
               end if
               this%max_collision = MAX(this%max_collision, ovflw_len)
            else
               tbl_entry%entry_value => newval
            end if
         END ASSOCIATE
      end if
      this%num_keys = this%num_keys + 1

   end subroutine hash_table_add_hash_key

end module hist_hash_table
