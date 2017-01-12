import StringIO

basic_types = [
("string256","character(256)"),
("string128","character(128)"),
("string80","character(80)"),
("string64","character(64)"),
("string32","character(32)"),
("string16","character(16)"),
("string8","character(8)"),
("real","real(4)"),
("double","real(8)"),
("int", "integer(4)"),
("long", "integer(8)"),
#Stupid fortran standard says you can't compare logicals with .eq.  gfortran follows this rule unless you use -fugly-logint
#Because that feature really held back C, didn't it?  There is even a snarky page on the gnu website explaining how it is really for the best.
#("logical", "logical"), 
]

def generate_list_functions(name,declaration):
	header =  """
type fl_{name}_list_item
	{declaration}, pointer :: value
	type(fl_{name}_list_item), pointer :: next
end type fl_{name}_list_item

type fl_{name}_list
	type(fl_{name}_list_item), pointer :: first, last, iter
	integer length
end type fl_{name}_list
""".format(name=name, declaration=declaration)

	body = """
subroutine fl_create_{name}_list(list)
	type(fl_{name}_list) :: list
	list%length=0
	nullify(list%first)
	nullify(list%last)
	nullify(list%iter)
end subroutine fl_create_{name}_list

subroutine fl_append_{name}_list(list,value)
	type(fl_{name}_list) :: list
	{declaration} :: value
	
	type(fl_{name}_list_item), pointer :: item
	allocate(item)
	allocate(item%value)
	item%value=value
	nullify(item%next)
	if (.not. associated(list%first)) then
		list%first => item
		list%last => item
		list%length=1
	else
		list%last%next => item
		list%last => item
		list%length=list%length+1
	endif
	call fl_reset_{name}_iterator_list(list)
end subroutine fl_append_{name}_list

subroutine fl_destroy_{name}_list(list)
	type(fl_{name}_list) :: list
	type(fl_{name}_list_item), pointer :: item, next

	if (.not. associated(list%first)) return
	item=>list%first
	do
		nullify(next)
		if (associated(item%next)) next=>item%next
		deallocate(item%value)
		deallocate(item)
		if (.not. associated(next)) exit
		item=>next
	enddo
	nullify(list%first)
	nullify(list%last)
	nullify(list%iter)
	list%length=0
end subroutine fl_destroy_{name}_list

function fl_{name}_list_to_array(list) result(arr)
	type(fl_{name}_list) :: list
	{declaration}, pointer, dimension(:) :: arr
	type(fl_{name}_list_item), pointer :: item
	integer i
	allocate(arr(list%length))
	item=>list%first
	if (list%length>0) then
		do i=1,list%length
			arr(i) = item%value
			if (associated(item%next)) item=>item%next
		enddo
	endif
end function fl_{name}_list_to_array

function fl_member_{name}_list(list,value) result(B)
	logical :: B
	type(fl_{name}_list) :: list
	type(fl_{name}_list_item), pointer :: item
	{declaration} :: value

	B=.false.
	if (associated(list%first)) then
		item => list%first
		do
			if (item%value == value) then
				B=.true.
				exit
			endif
			if (.not. associated(item%next)) exit
			item=>item%next
		enddo
	endif
end function fl_member_{name}_list

function fl_pop_last_{name}_list(list) result(pop)
	type(fl_{name}_list) :: list
	type(fl_{name}_list_item), pointer :: item	
	{declaration} :: pop
	integer i

	if (list%length==0) stop 'Requested last item from empty fl_{name}_list'
	
	pop=list%last%value
	deallocate(list%last%value)
	deallocate(list%last)
	list%length=list%length-1
	if (list%length==0) then
		nullify(list%first)
		nullify(list%last)
	elseif (list%length==1) then
		list%last=>list%first
	else
		item=>list%first
		do i=1,list%length-1
			item=>item%next
		enddo
		nullify(item%next)
	endif
	call fl_reset_{name}_iterator_list(list)	
end function fl_pop_last_{name}_list

function fl_pop_first_{name}_list(list) result(pop)
	type(fl_{name}_list) :: list
	type(fl_{name}_list_item), pointer :: item	
	{declaration} :: pop
	if (list%length==0) stop 'Requested first item from empty fl_{name}_list'
	item => list%first
	nullify(list%first)
	if (associated(item%next)) list%first=>item%next
	pop=item%value
	deallocate(item%value)
	deallocate(item)
	list%length=list%length-1
	call fl_reset_{name}_iterator_list(list)	
end function fl_pop_first_{name}_list

function fl_get_{name}_list(list,n) result(x)
	type(fl_{name}_list) :: list
	integer n,i
	type(fl_{name}_list_item), pointer :: item	
	{declaration} :: x
	if (n>list%length .or. n<1) stop 'Requested out-of-range fl_{name}_list item'
	item=>list%first
	do i=1,n-1
		item=>item%next
	enddo
	x = item%value
end function fl_get_{name}_list

function fl_iterate_{name}_list(list) result(x)
	type(fl_{name}_list) :: list
	{declaration} :: x
	if (list%length==0) stop 'Tried to iterate through empty fl_{name}_list'
	if (.not. associated(list%iter)) list%iter => list%first
	x=list%iter%value
	if (.not. associated(list%iter%next)) then
		nullify(list%iter)
	else
		list%iter => list%iter%next
	endif
end function fl_iterate_{name}_list

subroutine fl_insert_{name}_list(list,value,after)
	type(fl_{name}_list) :: list
	{declaration} :: value
	type(fl_{name}_list_item), pointer :: item, previous,next
	integer i
	integer after
	if (after>list%length .or. after<1) stop 'Tried to insert item in invalid position in fl_{name}_list'
	item=>list%first
	if (after==list%length) then
		item=>list%last
	else
		do i=1,after-1
			item=>item%next
		enddo
	endif
	previous=>item
	nullify(next)
	if (associated(item%next)) next=>item%next
	nullify(item)
	allocate(item)
	allocate(item%value)
	item%value=value
	previous%next=>item
	nullify(item%next)
	if (associated(item%next)) item%next=>next
	list%length = list%length + 1
	call fl_reset_{name}_iterator_list(list)	
end subroutine fl_insert_{name}_list

subroutine fl_reset_{name}_iterator_list(list)
	type(fl_{name}_list) :: list
	nullify(list%iter)
end subroutine fl_reset_{name}_iterator_list


""".format(name=name, declaration=declaration)
	return header,body

	
def generate_modules(types,module_name="fl_lists",uses=None):
	output=StringIO.StringIO()
	headers = []
	bodies = []
	for name,decl in types:
		header,body = generate_list_functions(name,decl)
		headers.append(header)
		bodies.append(body)
	output.write("module {0}\n".format(module_name))
	if uses:
		for use in uses:
			output.write("use {0}\n".format(use))
	output.write("implicit none\n")
	for header in headers:
		output.write("{0}\n".format(header))
	output.write("contains\n")
	for body in bodies:
		output.write("{0}\n".format(body))
	output.write("end module {0}\n".format(module_name))
	output.seek(0)
	return output.read()
		
	
if __name__=="__main__":
	print "!This module was auto-generated from build_list_module.py."
	print "!Edits to it will be lost."
	print
	print generate_modules(basic_types)