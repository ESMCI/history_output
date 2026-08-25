module hist_hashable

   implicit none
   private

   ! Public interfaces
   public :: new_hashable_char
   public :: new_hashable_int

   type, abstract, public :: hist_hashable_t
      ! The hashable type is a base type that contains a hash key.
   contains
      procedure(hist_hashable_get_key), deferred :: key
   end type hist_hashable_t

   type, public, extends(hist_hashable_t) :: hist_hashable_char_t
      character(len=:), private, allocatable :: name
   contains
      procedure :: key => hist_hashable_char_get_key
   end type hist_hashable_char_t

   type, public, extends(hist_hashable_t) :: hist_hashable_int_t
      integer, private :: value
   contains
      procedure :: key => hist_hashable_int_get_key
      procedure :: val => hist_hashable_int_get_val
   end type hist_hashable_int_t

   ! Abstract interface for key procedure of hist_hashable_t class
   abstract interface
      function hist_hashable_get_key(hashable) result(key)
         import :: hist_hashable_t
         implicit none
         class(hist_hashable_t), intent(in) :: hashable
         character(len=:), allocatable      :: key
      end function hist_hashable_get_key
   end interface

contains

   !#######################################################################

   subroutine new_hashable_char(name_in, new_obj)
      character(len=*), intent(in)                       :: name_in
      type(hist_hashable_char_t), pointer, intent(inout) :: new_obj

      if (associated(new_obj)) then
         deallocate(new_obj)
      end if
      allocate(new_obj)
      new_obj%name = name_in
   end subroutine new_hashable_char

   !#######################################################################

   function hist_hashable_char_get_key(hashable) result(key)
      ! Return the hashable char class key (name)
      class(hist_hashable_char_t), intent(in) :: hashable
      character(len=:), allocatable           :: key

      key = hashable%name
   end function hist_hashable_char_get_key

   !#######################################################################

   subroutine new_hashable_int(val_in, new_obj)
      integer,                            intent(in)    :: val_in
      type(hist_hashable_int_t), pointer, intent(inout) :: new_obj

      if (associated(new_obj)) then
         deallocate(new_obj)
      end if
      allocate(new_obj)
      new_obj%value = val_in
   end subroutine new_hashable_int

   !#######################################################################

   function hist_hashable_int_get_key(hashable) result(key)
      ! Return the hashable int class key (value ==> string)
      class(hist_hashable_int_t), intent(in) :: hashable
      character(len=:), allocatable          :: key

      character(len=32) :: key_str

      write(key_str, '(i0)') hashable%val()
      key = trim(key_str)
   end function hist_hashable_int_get_key

   !#######################################################################

   function hist_hashable_int_get_val(hashable) result(val)
      ! Return the hashable int class value
      class(hist_hashable_int_t), intent(in) :: hashable
      integer                                :: val

      val = hashable%value
   end function hist_hashable_int_get_val

end module hist_hashable
