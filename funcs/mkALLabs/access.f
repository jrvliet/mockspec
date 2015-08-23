	INTEGER FUNCTION access(file,type)

	character*(*) file, type
	open(99,file=file,status='old',iostat=ierr)

        access = ierr
        close(99)
	return
	end
