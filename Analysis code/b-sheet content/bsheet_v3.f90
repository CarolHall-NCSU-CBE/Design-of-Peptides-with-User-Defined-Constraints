!  bsheet_content.f90 
!

!****************************************************************************
!
!  PROGRAM: bsheet_content
!
!  Author:  Haoyu Wang.
!
! Need the "input.txt", "results" folder to work
!
!****************************************************************************

    module my_functions
    implicit none

    contains
    
    subroutine read_input(input_file, pep1, pep2, nc1, nc2, nb1, nb2)
    character(len=*) :: input_file
    character(len=256) :: pep1, pep2
    integer, intent(out) :: nc1, nc2, nb1, nb2

    character(len=256) :: line
    character(len=256) :: key, value, substr
    character(len=1) :: char
    integer :: i, ios, g_count, eqpos

    character(len=256), dimension(8) :: lines

    ! Read first 8 lines
    open(unit=10, file=input_file, status='old', action='read', iostat=ios)
    if (ios /= 0) then
        print *, 'Error opening file:', trim(input_file)
        stop
    end if

    do i = 1, 8
        read(10,'(A)', iostat=ios) line
        if (ios /= 0) exit
        lines(i) = adjustl(trim(line))
    end do
    close(10)

    ! read in pep1
    if (index(lines(2), '=') > 0) then
        call parse_line(lines(2), pep1)
    end if
    ! read in pep2
    if (index(lines(6), '=') > 0) then
        call parse_line(lines(6), pep2)
    end if

    ! read number of chain (nc1) for pep1
    eqpos = index(lines(4), '=')
    if (eqpos > 0) then
        substr = adjustl(trim(lines(4)(eqpos + 1:)))
        read(substr, *) nc1
    end if
    ! read number of chain (nc2) for pep2
    eqpos = index(lines(8), '=')
    if (eqpos > 0) then
        substr = adjustl(trim(lines(8)(eqpos + 1:)))
        read(substr, *) nc2
    end if

    ! Count beads: non-"G" = 4, "G" = 3
    g_count = 0
    do i = 1, len_trim(pep1)
        char = pep1(i:i)
        if (char == 'G') g_count = g_count + 1
    end do
    nb1 = 4 * len_trim(pep1) - g_count

    g_count = 0
    do i = 1, len_trim(pep2)
        char = pep2(i:i)
        if (char == 'G') g_count = g_count + 1
    end do
    nb2 = 4 * len_trim(pep2) - g_count

    end subroutine read_input

    subroutine parse_line(line, outputs)
        implicit none
        character(len=*), intent(in) :: line
        character(len=*), intent(out) :: outputs
        integer :: eq_pos

        eq_pos = index(line, '=')
        if (eq_pos > 0) then
            outputs = adjustl(trim(line(eq_pos + 1:)))
        else
            outputs = ''
        end if
    end subroutine parse_line

    subroutine read_pdb(filename, atom_name, aa_name, resid, xyz, chain_id, natoms, nchain)
    implicit none

    character(len=*), intent(in)                 :: filename
    character(len=4), allocatable, intent(out)   :: atom_name(:)
    character(len=3), allocatable, intent(out)   :: aa_name(:)
    integer, allocatable, intent(out)            :: resid(:), chain_id(:)
    real(8), allocatable, intent(out)            :: xyz(:,:)
    integer, intent(out)                         :: natoms, nchain

    character(len=256)                           :: line
    integer                                      :: ios, last_resid

    natoms = 0
    open(unit=11, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) then
        write(*,'(2A)') 'Error opening PDB file: ', trim(filename)
        stop
    endif

    do
        read(11,'(A)',iostat=ios) line
        if (ios /= 0) exit
        if (line(1:6) == 'ATOM  ' .or. line(1:6) == 'HETATM') natoms = natoms + 1
    enddo

    if (natoms == 0) then
        write(*,'(2A)') 'No atoms were found in PDB file: ', trim(filename)
        close(11)
        stop
    endif

    allocate(atom_name(natoms), aa_name(natoms), resid(natoms), chain_id(natoms), xyz(3,natoms))

    rewind(11)
    natoms = 0
    nchain = 0
    last_resid = 0

    do
        read(11,'(A)',iostat=ios) line
        if (ios /= 0) exit
        if (line(1:6) /= 'ATOM  ' .and. line(1:6) /= 'HETATM') cycle

        natoms = natoms + 1
        atom_name(natoms) = adjustl(line(13:16))
        aa_name(natoms) = adjustl(line(18:20))
        read(line(23:26),*,iostat=ios) resid(natoms)
        if (ios /= 0) then
            write(*,'(A,I0,A)') 'Cannot read residue number from atom ', natoms, ' in the PDB file.'
            close(11)
            stop
        endif
        read(line(31:38),*,iostat=ios) xyz(1,natoms)
        if (ios == 0) read(line(39:46),*,iostat=ios) xyz(2,natoms)
        if (ios == 0) read(line(47:54),*,iostat=ios) xyz(3,natoms)
        if (ios /= 0) then
            write(*,'(A,I0,A)') 'Cannot read coordinates from atom ', natoms, ' in the PDB file.'
            close(11)
            stop
        endif

        if (natoms == 1 .or. resid(natoms) < last_resid .or. &
            (resid(natoms) == 1 .and. trim(atom_name(natoms)) == 'N')) then
            nchain = nchain + 1
        endif
        chain_id(natoms) = nchain
        last_resid = resid(natoms)
    enddo

    close(11)
    end subroutine read_pdb


    subroutine analyze_pdb(pep, nc, nb, nres, numseq, pdb_file)
    implicit none

    character(len=256), intent(out)               :: pep(10)
    integer, intent(out)                          :: nc(10), nb(10), nres(10), numseq
    character(len=100), intent(out)               :: pdb_file

    character(len=4), allocatable                 :: atom_name(:)
    character(len=3), allocatable                 :: aa_name(:)
    character(len=256), allocatable               :: chain_sequence(:)
    character(len=1)                              :: one_letter
    integer, allocatable                          :: resid(:), chain_id(:), last_ca_resid(:)
    real(8), allocatable                          :: xyz(:,:)
    integer                                       :: natoms, nchain, irun, i, j, k, ngly
    logical                                       :: pdb_exists, sequence_found

    pdb_exists = .false.
    do irun = 0, 9999
        write(pdb_file,'(A,I4.4,A)') 'results/run', irun, '.pdb'
        inquire(file=trim(pdb_file), exist=pdb_exists)
        if (pdb_exists) exit
    enddo
    if (.not. pdb_exists) then
        write(*,'(A)') 'No PDB file matching results/run0000.pdb through results/run9999.pdb was found.'
        stop
    endif

    call read_pdb(pdb_file, atom_name, aa_name, resid, xyz, chain_id, natoms, nchain)

    allocate(chain_sequence(nchain), last_ca_resid(nchain))
    chain_sequence = ''
    last_ca_resid = -1

    do i = 1, natoms
        if (trim(atom_name(i)) /= 'CA') cycle
        if (resid(i) == last_ca_resid(chain_id(i))) cycle

        select case (trim(aa_name(i)))
        case ('ALA'); one_letter = 'A'
        case ('ARG'); one_letter = 'R'
        case ('ASN'); one_letter = 'N'
        case ('ASP'); one_letter = 'D'
        case ('CYS'); one_letter = 'C'
        case ('GLN'); one_letter = 'Q'
        case ('GLU'); one_letter = 'E'
        case ('GLY'); one_letter = 'G'
        case ('HIS'); one_letter = 'H'
        case ('ILE'); one_letter = 'I'
        case ('LEU'); one_letter = 'L'
        case ('LYS'); one_letter = 'K'
        case ('MET'); one_letter = 'M'
        case ('PHE'); one_letter = 'F'
        case ('PRO'); one_letter = 'P'
        case ('SER'); one_letter = 'S'
        case ('THR'); one_letter = 'T'
        case ('TRP'); one_letter = 'W'
        case ('TYR'); one_letter = 'Y'
        case ('VAL'); one_letter = 'V'
        case default
            write(*,'(3A,I0,A)') 'Unknown residue ', trim(aa_name(i)), ' in chain ', chain_id(i), '.'
            stop
        end select

        j = len_trim(chain_sequence(chain_id(i))) + 1
        if (j > len(chain_sequence(chain_id(i)))) then
            write(*,'(A,I0,A)') 'Peptide chain ', chain_id(i), ' is longer than the sequence array.'
            stop
        endif
        chain_sequence(chain_id(i))(j:j) = one_letter
        last_ca_resid(chain_id(i)) = resid(i)
    enddo

    pep = ''
    nc = 0
    nb = 0
    nres = 0
    numseq = 0

    do i = 1, nchain
        if (len_trim(chain_sequence(i)) == 0) cycle
        sequence_found = .false.
        do j = 1, numseq
            if (trim(chain_sequence(i)) == trim(pep(j))) then
                nc(j) = nc(j) + 1
                sequence_found = .true.
                exit
            endif
        enddo

        if (.not. sequence_found) then
            if (numseq == 10) then
                write(*,'(A)') 'More than 10 unique peptide sequences were found in the PDB file.'
                stop
            endif
            numseq = numseq + 1
            pep(numseq) = trim(chain_sequence(i))
            nc(numseq) = 1
        endif
    enddo

    if (numseq == 0) then
        write(*,'(A)') 'No peptide sequences could be constructed from the CA records in the PDB file.'
        stop
    endif

    do i = 1, numseq
        nres(i) = len_trim(pep(i))
        ngly = 0
        do k = 1, nres(i)
            if (pep(i)(k:k) == 'G') ngly = ngly + 1
        enddo
        nb(i) = 4*nres(i) - ngly
    enddo

    write(*,'(2A)') 'PDB file used: ', trim(pdb_file)
    write(*,'(A,I0)') 'Number of unique peptides: ', numseq
    do i = 1, numseq
        write(*,'(A,I0,2A)') 'Pep', i, ': ', trim(pep(i))
        write(*,'(A,I0)') '  Number of chains : ', nc(i)
        write(*,'(A,I0)') '  Number of beads : ', nb(i)
        write(*,'(A)') " "
    enddo

    deallocate(atom_name, aa_name, resid, xyz, chain_id, chain_sequence, last_ca_resid)
    end subroutine analyze_pdb
    
    
    subroutine read_energy(energy_file, ttotal, colltotal, i, max_run, anneal_steps, &
        t_series, coll_series, npts)
    implicit none
    character(len=*)                  :: energy_file
    real(8)                           :: ttotal, timetemp, factor
    real(8)                           :: dummy1, dummy2, dummy3, dummy4
    real(8)                           :: dummy5, dummy6, dummy7, dummy8, dummy9
    integer, intent(in)               :: i, max_run
    integer, intent(inout)            :: anneal_steps
    integer                           :: npts, unit, ios, count, j, anneal_line
    integer(8), intent(inout)         :: colltotal
    integer(8)                        :: anneal_collisions, dummy0
    real(8), allocatable              :: t_series(:), t(:)
    integer(8), allocatable           :: coll_series(:), coll(:)
    character(len=100)                :: line, temp_file
    logical                           :: new_anneal_style, file_exists
    logical                           :: temperature_found

    ! Determine the annealing input style and number of steps once.
    inquire(file='inputs/annealtemps', exist=new_anneal_style)

    if (anneal_steps < 0) then
        anneal_steps = 0

        if (new_anneal_style) then
            open(10, file='inputs/annealtemps', status='old', action='read', iostat=ios)

            do
                read(10,'(A)',iostat=ios) line
                if (ios /= 0) exit
                line = adjustl(line)
                if (len_trim(line) == 0) cycle
                if (line(1:1) == '#') cycle

                read(line,*,iostat=ios) timetemp, anneal_collisions
                anneal_steps = anneal_steps + 1
            enddo
            close(10)
            write(*,'(A)') 'Annealing input file: inputs/annealtemps'
        else
            do j = 1, 9999
                write(temp_file,'(A,I0)') 'inputs/annealtemp_', j
                inquire(file=trim(temp_file), exist=file_exists)
                if (.not. file_exists) exit
                anneal_steps = j
            enddo

            if (anneal_steps > 0) then
                write(*,'(A)') 'Annealing input files: inputs/annealtemp_1, 2, 3...'
            else
                write(*,'(A)') 'No annealing-temperature input was found.'
            endif
        endif

        write(*,'(A,I0)') 'Number of annealing steps: ', anneal_steps

        if (max_run < anneal_steps) then
            anneal_steps = max_run
        endif
    endif

    ! Determine the reduced temperature for this run.
    if (i == 0) then
        timetemp = 0.5d0
    elseif (i <= anneal_steps) then
        if (new_anneal_style) then
            open(10, file='inputs/annealtemps', status='old', action='read', iostat=ios)
            if (ios /= 0) then
                write(*,'(A)') 'Cannot open the new-style annealing file: inputs/annealtemps'
                stop
            endif

            anneal_line = 0
            temperature_found = .false.
            do
                read(10,'(A)',iostat=ios) line
                if (ios /= 0) exit
                line = adjustl(line)
                if (len_trim(line) == 0) cycle
                if (line(1:1) == '#') cycle

                anneal_line = anneal_line + 1
                if (anneal_line == i) then
                    read(line,*,iostat=ios) timetemp, anneal_collisions
                    temperature_found = .true.
                    exit
                endif
            enddo
            close(10)

        else
            write(temp_file, '(A,I0)') 'inputs/annealtemp_', i
            open(10, file=trim(temp_file), status='old', action='read', iostat=ios)
            read(10,*,iostat=ios) timetemp
            close(10)
        endif
    else
        temp_file = 'inputs/simtemp'
        open(10, file=trim(temp_file), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(*,'(A)') 'Cannot open simulation-temperature file: inputs/simtemp'
            stop
        endif
        read(10,*,iostat=ios) timetemp
        close(10)
    end if

    ! 2. Read energy file, we only read the reduce time
    open(10, file=trim(energy_file), status='old', action='read')
    allocate(t(10000), coll(10000))  ! Max lines (adjust as needed)
    count = 0
    
    do
        read(10, *, iostat=ios) dummy0, t(count+1), dummy1, dummy2, dummy3, dummy4, &
            dummy5, dummy6, dummy7, dummy8, dummy9
        if (ios /= 0) exit
        count = count + 1
        coll(count) = dummy0
    end do
    close(10)

    ! count is length of time array, npts is length the the array remove the first one (which is t=0)
    if (count <= 1) then
        allocate(t_series(0), coll_series(0))
        npts = 0
        return
    end if

    allocate(t_series(count - 1), coll_series(count - 1))
    do j = 2, count
        t_series(j - 1) = ttotal + t(j) * 0.96d0 * 0.001d0 * 3.3d0 / sqrt(timetemp * 12.0d0)
        coll_series(j - 1) = colltotal + coll(j)
    end do
    ttotal = t_series(count - 1)
    colltotal = coll_series(count - 1)
    npts = count - 1

    end subroutine read_energy
    
    subroutine read_bptnr_multi(filename, nc1, nc2, nb1, nb2, nframes, collision_nums, bptnrs, hbmats)
    implicit none
    character(len=*) :: filename
    integer                                 :: nc1, nc2, nb1, nb2, nframes
    integer                                 :: nop1, nop2, noptotal, nchain
    integer                                 :: unit, ios, i, j, frame    
    integer                                 :: bead, partner, chain_i, chain_j
    
    integer, allocatable                    :: bptnrs(:,:,:)
    integer, allocatable                    :: hbmats(:,:,:)
    integer, allocatable                    :: hb(:)    
    integer(8), allocatable                 :: collision_nums(:)
    integer(8)                              :: coll
    
    ! Initialize
    nop1 = nc1 * nb1
    nop2 = nc2 * nb2
    noptotal = nop1 + nop2
    nchain = nc1 + nc2

    allocate(hb(noptotal))
    allocate(collision_nums(0))
    allocate(bptnrs(noptotal, 2, 0))
    allocate(hbmats(nchain, nchain, 0))

    frame = 0
    open(unit=15, file=filename, form='unformatted', status='old', action='read', iostat=ios)
    if (ios /= 0) then
        print *, 'Error opening file: ', trim(filename)
        stop
    end if

    do
        read(15, iostat=ios) coll, hb
        if (ios /= 0) exit

        frame = frame + 1

        ! Expand arrays
        call expand_arrays(frame, noptotal, nchain, collision_nums, bptnrs, hbmats)

        collision_nums(frame) = coll

        do i = 1, noptotal
            bptnrs(i, 1, frame) = i
            bptnrs(i, 2, frame) = hb(i)
        end do

        hbmats(:, :, frame) = 0
        do i = 1, noptotal
            partner = hb(i)
            if (partner /= 0) then
                if (i <= nop1) then
                    chain_i = (i - 1) / nb1 + 1
                else
                    chain_i = (i - nop1 - 1) / nb2 + nc1 + 1
                end if

                if (partner <= nop1) then
                    chain_j = (partner - 1) / nb1 + 1
                else
                    chain_j = (partner - nop1 - 1) / nb2 + nc1 + 1
                end if

                hbmats(chain_i, chain_j, frame) = hbmats(chain_i, chain_j, frame) + 1
            end if
        end do
    end do
    close(15)

    nframes = frame
    
    deallocate(hb)
    end subroutine read_bptnr_multi
    
    
    subroutine expand_arrays(frame, noptotal, nchain, collision_nums, bptnrs, hbmats)
    implicit none
    integer, intent(in) :: frame, noptotal, nchain
    integer(8), allocatable, intent(inout) :: collision_nums(:)
    integer, allocatable, intent(inout) :: bptnrs(:,:,:)
    integer, allocatable, intent(inout) :: hbmats(:,:,:)

    integer(8), allocatable :: temp_col(:)
    integer, allocatable :: temp_bptnr(:,:,:)
    integer, allocatable :: temp_hbm(:,:,:)

    allocate(temp_col(frame))
    if (frame > 1) temp_col(1:frame-1) = collision_nums
    deallocate(collision_nums)
    collision_nums = temp_col

    allocate(temp_bptnr(noptotal, 2, frame))
    if (frame > 1) temp_bptnr(:,:,1:frame-1) = bptnrs
    deallocate(bptnrs)
    bptnrs = temp_bptnr

    allocate(temp_hbm(nchain, nchain, frame))
    if (frame > 1) temp_hbm(:,:,1:frame-1) = hbmats
    deallocate(hbmats)
    hbmats = temp_hbm
    end subroutine expand_arrays

    
    
    subroutine bsheet_content(bptnr, hbmatrix, nc1, nc2, nb1, nb2, nres1, nres2, anti, para, turn, undetermined, content)
    implicit none
    ! Inputs
    integer, intent(in) :: nc1, nc2, nb1, nb2, nres1, nres2
    integer, intent(in) :: bptnr(:,:), hbmatrix(:,:)

    ! Outputs
    integer, intent(out) :: anti, para, turn, undetermined
    real(8), intent(out) :: content

    ! Locals
    integer :: nchain, total_residues
    integer, allocatable                    :: bsheet_assign(:)
    integer                                 :: i, j, m, ii, jj, jj2, ii2
    integer                                 :: nbi, nbii, nresi, nbj, nbjj, nresj
    integer, allocatable                    :: nh_i(:), nh_j(:), co_i(:), co_j(:)  
    integer, allocatable                    :: resNHi(:), resCOi(:), resNHj(:), resCOj(:) 
    integer                                 :: nresii, nresjj, nbond_ij, nbond_ji
    integer                                 :: NH_i_pos, CO_i_pos, CO_j_pos, NH_j_pos
    integer                                 :: k, x, counter, diff, nh_basg, co_basg
    real(8)                                 :: a, b
    logical                                 :: found

    anti = 0; para = 0; turn = 0; undetermined = 0
    nchain = nc1 + nc2

    total_residues = nc1 * nres1 + nc2 * nres2
    allocate(bsheet_assign(total_residues))
    bsheet_assign = 0

    do i = 1, nchain
        do j = i, nchain
            nbond_ij = 0; nbond_ji = 0                ! nbond_i: number of bonds (NH_i to CO_j), and same ......
            if (hbmatrix(i, j) > 0) then
                ! Determine nb and nres for chains i and j
                if (i <= nc1) then                      ! if current i is pep1
                    nbi = nb1                           ! number of beads in pep1
                    nresi = nres1                       ! number of residues in pep1
                    nbii = (i - 1) * nb1 + 1            ! chain i bead starting index
                    nresii = (i - 1) * nres1            ! chain i residue starting index
                else
                    nbi = nb2
                    nresi = nres2
                    nbii = (i - nc1 - 1) * nb2 + nc1 * nb1  + 1
                    nresii = (i - nc1 - 1) * nres2 + nc1 * nres1
                endif

                if (j <= nc1) then
                    nbj = nb1
                    nresj = nres1
                    nbjj = (j - 1) * nb1 + 1
                    nresjj = (j - 1) * nres1
                else
                    nbj = nb2
                    nresj = nres2
                    nbjj = (j - nc1 - 1) * nb2 + nb1 * nc1 + 1
                    nresjj = (j - nc1 - 1) * nres2 + nres1 * nc1
                endif

                
                allocate(nh_i(nresi))
                allocate(nh_j(nresj))
                allocate(co_i(nresi))
                allocate(co_j(nresj))
                nh_i = 0; co_j = 0; nh_j = 0; co_i = 0

                do m = 1, nresi
                    CO_j_pos = bptnr(nbii + nresi + m - 1, 2)
                    
                    jj  = int(real(CO_j_pos/ nb1 + 1))
                    jj2 = int(real(CO_j_pos  - nc1*nb1)/ nb2 + 1 + nc1) 
                    
                    if (( jj == j .or. jj2 == j) .and. (CO_j_pos /=0 )) then
                        nh_i(m) = nbii + nresi + m - 1             ! index of m_th NH bead in chain i
                        co_j(m) = CO_j_pos                         ! index of m_th CO bead in chain j
                        nbond_ij = nbond_ij + 1
                    endif
                enddo

                     
                do m = 1, nresj
                    CO_i_pos = bptnr(nbjj + nresj + m - 1, 2)
                    ii  = int(real(CO_i_pos/ nb1 + 1))
                    ii2 = int(real(CO_i_pos  - nc1*nb1)/ nb2 + 1 + nc1)    
                    if (( ii == i .or. ii2 == i) .and. (CO_i_pos /=0 )) then
                        nh_j(m) = nbjj + nresj + m - 1
                        co_i(m) = CO_i_pos
                        nbond_ji = nbond_ji + 1

                        
                    endif
                enddo

                nh_i = pack(nh_i, nh_i /= 0)
                nh_j = pack(nh_j, nh_j /= 0)
                co_i = pack(co_i, co_i /= 0)
                co_j = pack(co_j, co_j /= 0)
                
                
                ! Determine sheet direction from CO changes
                a = 0.0
                b = 0.0
                if (nbond_ij >= 2) then
                    do m = 2, nbond_ij
                        if (co_j(m) /= 0 .and. co_j(m - 1) /= 0) then
                            a = a + (co_j(m) - co_j(m - 1))
                        endif
                    enddo
                endif

                if (nbond_ji >= 2) then
                    do m = 2, nbond_ji
                        if (co_i(m) /= 0 .and. co_i(m - 1) /= 0) then
                            b = b + (co_i(m) - co_i(m - 1))
                        endif
                    enddo
                endif

                if (i <= nc1) then                  
                    resNHi = nh_i - nresi - (i - 1) * nbi                   ! calculate NH is which residue in chain i, array calculation
                    resCOi = co_i - 2*nresi - (i - 1) * nbi                 ! calculate CO is which residue in chain i
                else
                    resNHi = nh_i - nresi - (i - nc1 - 1) * nbi - nc1 * nb1
                    resCOi = co_i - 2*nresi - (i - nc1 - 1) * nbi - nc1 * nb1
                endif

                if (j <= nc1) then
                    resCOj = co_j - 2*nresj - (j - 1) * nbj
                    resNHj = nh_j - nresj - (j - 1) * nbj
                else
                    resCOj = co_j - 2*nresj - (j - nc1 - 1) * nbj - nc1 * nb1
                    resNHj = nh_j - nresj - (j - nc1 - 1) * nbj - nc1 * nb1
                endif
                
!!!!!!!!!!!!!!!!!!!!! LOOP 1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                do m = 1, size(resNHi)
                    NH_i_pos = resNHi(m)
                    CO_j_pos = resCOj(m)

                    ! Check if CO_j_pos exists in resNHj
                    found = .false.
                    do k = 1, size(resNHj)
                        if (resNHj(k) == CO_j_pos) then
                            found = .true.
                            CO_i_pos = resCOi(k)

                            if (CO_i_pos == NH_i_pos) then                      ! antiparallel
                                if (bsheet_assign(nresii + NH_i_pos) <= 1) then
                                    bsheet_assign(nresii + NH_i_pos) = 2
                                else if (bsheet_assign(nresii + NH_i_pos) == 3) then
                                    bsheet_assign(nresii + NH_i_pos) = 5
                                endif

                                if (bsheet_assign(nresjj + CO_j_pos) <= 1) then
                                    bsheet_assign(nresjj + CO_j_pos) = 2
                                else if (bsheet_assign(nresjj + CO_j_pos) == 3) then
                                    bsheet_assign(nresjj + CO_j_pos) = 5
                                endif

                            else if (CO_i_pos == NH_i_pos - 2) then             ! parallel
                                if (bsheet_assign(nresii + NH_i_pos) <= 1) then
                                    bsheet_assign(nresii + NH_i_pos) = 3
                                else if (bsheet_assign(nresii + NH_i_pos) == 2) then
                                    bsheet_assign(nresii + NH_i_pos) = 5
                                endif

                                if(NH_i_pos - 1 > 0) then
                                    if (bsheet_assign(nresii + NH_i_pos - 1) <= 1) then
                                        bsheet_assign(nresii + NH_i_pos - 1) = 3
                                    else if (bsheet_assign(nresii + NH_i_pos - 1) == 2) then
                                        bsheet_assign(nresii + NH_i_pos - 1) = 5
                                    endif
                                endif

                                if (bsheet_assign(nresjj + CO_j_pos) <= 1) then
                                    bsheet_assign(nresjj + CO_j_pos) = 3
                                else if (bsheet_assign(nresjj + CO_j_pos) == 2) then
                                    bsheet_assign(nresjj + CO_j_pos) = 5
                                endif

                                if (bsheet_assign(nresii + CO_i_pos) <= 1) then
                                    bsheet_assign(nresii + CO_i_pos) = 3
                                else if (bsheet_assign(nresii + CO_i_pos) == 2) then
                                    bsheet_assign(nresii + CO_i_pos) = 5
                                end if
                            end if
                            
                            exit
                        end if
                    end do
                
                    if (.not. found) then
                        ! CO_j_pos + 2 in resNHj
                        found = .false.
                        do k = 1, size(resNHj)
                            if (resNHj(k) == CO_j_pos + 2) then
                                found = .true.
                                CO_i_pos = resCOi(k)

                                if (CO_i_pos == NH_i_pos) then             ! parallel
                                    if (bsheet_assign(nresii + NH_i_pos) <= 1) then
                                        bsheet_assign(nresii + NH_i_pos) = 3
                                    else if (bsheet_assign(nresii + NH_i_pos) == 2) then
                                        bsheet_assign(nresii + NH_i_pos) = 5
                                    end if
                
                                    if (bsheet_assign(nresjj + CO_j_pos) <= 1) then
                                        bsheet_assign(nresjj + CO_j_pos) = 3
                                    else if (bsheet_assign(nresjj + CO_j_pos) == 2) then
                                        bsheet_assign(nresjj + CO_j_pos) = 5
                                    end if
                
                                    if (bsheet_assign(nresjj + CO_j_pos + 1) <= 1) then
                                        bsheet_assign(nresjj + CO_j_pos + 1) = 3
                                    else if (bsheet_assign(nresjj + CO_j_pos + 1) == 2) then
                                        bsheet_assign(nresjj + CO_j_pos + 1) = 5
                                    end if

                                    if (bsheet_assign(nresjj + CO_j_pos + 2) <= 1) then
                                        bsheet_assign(nresjj + CO_j_pos + 2) = 3
                                    else if (bsheet_assign(nresjj + CO_j_pos + 2) == 2) then
                                        bsheet_assign(nresjj + CO_j_pos + 2) = 5
                                    end if
                                end if
                                exit
                            end if
                        end do

                        if (.not. found) then
                            ! no pair — single H-bond
                            if (bsheet_assign(nresii + NH_i_pos) <= 1) then
                                bsheet_assign(nresii + NH_i_pos) = 1
                            end if
                            
                            if (bsheet_assign(nresjj + CO_j_pos) <= 1) then
                                bsheet_assign(nresjj + CO_j_pos) = 1
                            end if
                        end if
                    end if
                            
                    found = .false.
                    do k = 1, size(resCOi)
                        if (resCOi(k) == NH_i_pos - 2) then
                            found = .true.
                            NH_j_pos = resCOi(k)  ! same as CO_i_p2_pos
                            exit
                        end if
                    end do

                    if (found) then
                        if (NH_j_pos == CO_j_pos + 2) then
                            if (bsheet_assign(nresii + NH_i_pos) <= 1) then
                                bsheet_assign(nresii + NH_i_pos) = 2
                            else if (bsheet_assign(nresii + NH_i_pos) == 3) then
                                bsheet_assign(nresii + NH_i_pos) = 5
                            end if

                            if(NH_i_pos - 1 > 0) then
                                if (bsheet_assign(nresii + NH_i_pos - 1) <= 1) then
                                    bsheet_assign(nresii + NH_i_pos - 1) = 2
                                else if (bsheet_assign(nresii + NH_i_pos - 1) == 3) then
                                    bsheet_assign(nresii + NH_i_pos - 1) = 5
                                endif
                            endif

                            if (NH_i_pos - 2 > 0) then
                                if (bsheet_assign(nresii + NH_i_pos - 2) <= 1) then
                                    bsheet_assign(nresii + NH_i_pos - 2) = 2
                                else if (bsheet_assign(nresii + NH_i_pos - 2) == 3) then
                                    bsheet_assign(nresii + NH_i_pos - 2) = 5
                                endif
                            endif

                            if (bsheet_assign(nresjj + CO_j_pos) <= 1) then
                                bsheet_assign(nresjj + CO_j_pos) = 2
                            else if (bsheet_assign(nresjj + CO_j_pos) == 3) then
                                bsheet_assign(nresjj + CO_j_pos) = 5
                            endif

                            if (bsheet_assign(nresjj + CO_j_pos + 1) <= 1) then
                                bsheet_assign(nresjj + CO_j_pos + 1) = 2
                            else if (bsheet_assign(nresjj + CO_j_pos + 1) == 3) then
                                bsheet_assign(nresjj + CO_j_pos + 1) = 5
                            end if

                            if (bsheet_assign(nresjj + CO_j_pos + 2) <= 1) then
                                bsheet_assign(nresjj + CO_j_pos + 2) = 2
                            else if (bsheet_assign(nresjj + CO_j_pos + 2) == 3) then
                                bsheet_assign(nresjj + CO_j_pos + 2) = 5
                            end if

                        end if
                    end if     
                end do
!!!!!!!!!!!!!!!!!!!!! LOOP 2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                
                do m = 1, size(resNHj)
                    NH_j_pos = resNHj(m)
                    CO_i_pos = resCOi(m)

                    ! Check if CO_j_pos exists in resNHj
                    found = .false.
                    do k = 1, size(resNHi)
                        if (resNHi(k) == CO_i_pos) then
                            found = .true.
                            CO_j_pos = resCOj(k)

                            if (CO_j_pos == NH_j_pos) then                      ! antiparallel
                                if (bsheet_assign(nresjj + NH_j_pos) <= 1) then
                                    bsheet_assign(nresjj + NH_j_pos) = 2
                                else if (bsheet_assign(nresjj + NH_j_pos) == 3) then
                                    bsheet_assign(nresjj + NH_j_pos) = 5
                                end if

                                if (bsheet_assign(nresii + CO_i_pos) <= 1) then
                                    bsheet_assign(nresii + CO_i_pos) = 2
                                else if (bsheet_assign(nresii + CO_i_pos) == 3) then
                                    bsheet_assign(nresii + CO_i_pos) = 5
                                end if
!!!!
                            else if (CO_i_pos == NH_i_pos - 2) then             ! parallel
                                if (bsheet_assign(nresjj + NH_j_pos) <= 1) then
                                    bsheet_assign(nresjj + NH_j_pos) = 3
                                else if (bsheet_assign(nresjj + NH_j_pos) == 2) then
                                    bsheet_assign(nresjj + NH_j_pos) = 5
                                end if

                                if (NH_j_pos - 1 > 0) then
                                    if (bsheet_assign(nresjj + NH_j_pos - 1) <= 1) then
                                        bsheet_assign(nresjj + NH_j_pos - 1) = 3
                                    else if (bsheet_assign(nresjj + NH_j_pos - 1) == 2 ) then
                                        bsheet_assign(nresjj + NH_j_pos - 1) = 5
                                    end if
                                endif

                                if (bsheet_assign(nresii + CO_i_pos) <= 1) then
                                    bsheet_assign(nresii + CO_i_pos) = 3
                                else if (bsheet_assign(nresii + CO_i_pos) == 2) then
                                    bsheet_assign(nresii + CO_i_pos) = 5
                                end if

                                if (bsheet_assign(nresjj + CO_j_pos) <= 1) then
                                    bsheet_assign(nresjj + CO_j_pos) = 3
                                else if (bsheet_assign(nresjj + CO_j_pos) == 2) then
                                    bsheet_assign(nresjj + CO_j_pos) = 5
                                end if
                            end if
                            
                            exit
                        end if
                    end do
!!!!!                
                    if (.not. found) then
                        found = .false.
                        do k = 1, size(resNHi)
                            if (resNHi(k) == CO_i_pos + 2) then
                                found = .true.
                                CO_j_pos = resCOj(k)

                                if (CO_j_pos == NH_j_pos) then             ! parallel
                                    if (bsheet_assign(nresjj + NH_j_pos) <= 1) then
                                        bsheet_assign(nresjj + NH_j_pos) = 3
                                    else if (bsheet_assign(nresjj + NH_j_pos) == 2) then
                                        bsheet_assign(nresjj + NH_j_pos) = 5
                                    end if
                
                                    if (bsheet_assign(nresii + CO_i_pos) <= 1) then
                                        bsheet_assign(nresii + CO_i_pos) = 3
                                    else if (bsheet_assign(nresii + CO_i_pos) == 2) then
                                        bsheet_assign(nresii + CO_i_pos) = 5
                                    end if
                
                                    if (bsheet_assign(nresii + CO_i_pos + 1) <= 1) then
                                        bsheet_assign(nresii + CO_i_pos + 1) = 3
                                    else if (bsheet_assign(nresii + CO_i_pos + 1) == 2) then
                                        bsheet_assign(nresii + CO_i_pos + 1) = 5
                                    end if

                                    if (bsheet_assign(nresii + CO_i_pos + 2) <= 1) then
                                        bsheet_assign(nresii + CO_i_pos + 2) = 3
                                    else if (bsheet_assign(nresii + CO_i_pos + 2) == 2) then
                                        bsheet_assign(nresii + CO_i_pos + 2) = 5
                                    end if
                                end if
                                exit
                            end if
                        end do

                        if (.not. found) then
                            ! no pair — single H-bond
                            if (bsheet_assign(nresjj + NH_j_pos) <= 1) then
                                bsheet_assign(nresjj + NH_j_pos) = 1
                            end if
                            
                            if (bsheet_assign(nresii + CO_i_pos) <= 1) then
                                bsheet_assign(nresii + CO_i_pos) = 1
                            end if
                        end if
                    end if
                            
                    found = .false.
                    do k = 1, size(resCOj)
                        if (resCOj(k) == NH_j_pos - 2) then
                            found = .true.
                            NH_i_pos = resCOj(k)  ! same as CO_i_p2_pos
                            exit
                        end if
                    end do

                    if (found) then
                        if (NH_i_pos == CO_i_pos + 2) then
                            if (bsheet_assign(nresjj + NH_j_pos) <= 1) then
                                bsheet_assign(nresjj + NH_j_pos) = 2
                            else if (bsheet_assign(nresjj + NH_j_pos) == 3) then
                                bsheet_assign(nresjj + NH_j_pos) = 5
                            end if
                            
                            if(NH_j_pos - 1 > 0) then
                                if (bsheet_assign(nresjj + NH_j_pos - 1) <= 1) then
                                    bsheet_assign(nresjj + NH_j_pos - 1) = 2
                                else if (bsheet_assign(nresjj + NH_j_pos - 1) == 3) then
                                    bsheet_assign(nresjj + NH_j_pos - 1) = 5
                                endif
                            endif
                            
                            if (NH_j_pos - 2 > 0) then
                                if (bsheet_assign(nresjj + NH_j_pos - 2) <= 1) then
                                    bsheet_assign(nresjj + NH_j_pos - 2) = 2
                                else if (bsheet_assign(nresjj + NH_j_pos - 2) == 3) then
                                    bsheet_assign(nresjj + NH_j_pos - 2) = 5
                                endif
                            endif
                            
                            if (bsheet_assign(nresii + CO_i_pos) <= 1) then
                                bsheet_assign(nresii + CO_i_pos) = 2
                            else if (bsheet_assign(nresii + CO_i_pos) == 3) then
                                bsheet_assign(nresii + CO_i_pos) = 5
                            endif

                            if (bsheet_assign(nresii + CO_i_pos + 1) <= 1) then
                                bsheet_assign(nresii + CO_i_pos + 1) = 2
                            else if (bsheet_assign(nresii + CO_i_pos + 1) == 3) then
                                bsheet_assign(nresii + CO_i_pos + 1) = 5
                            end if

                            if (bsheet_assign(nresii + CO_i_pos + 2) <= 1) then
                                bsheet_assign(nresii + CO_i_pos + 2) = 2
                            else if (bsheet_assign(nresii + CO_i_pos + 2) == 3) then
                                bsheet_assign(nresii + CO_i_pos + 2) = 5
                            end if

                        end if
                    end if     
                end do                
                
                
                if ((a>0 .and. b>=0) .or. (a>=0 .and. b>0)) then        ! the two chains are parallel
                    do x=1, nresi                                       ! two chains are equal length, otherwise, need two judgement
                        if ((bsheet_assign(nresii + x) == 3) .or. (bsheet_assign(nresii + x) == 5)) para = para + 1
                    enddo
                    
                    do x=1, nresj
                        if ((bsheet_assign(nresjj + x) == 3) .or. (bsheet_assign(nresjj + x) == 5)) para = para + 1
                    enddo
                else if (( a<0 .and. b<=0) .or. (a<=0 .and. b<0)) then
                    do x=1, nresj
                        if ((bsheet_assign(nresjj + x) == 2) .or. (bsheet_assign(nresjj + x) == 5)) anti = anti + 1
                    enddo
                    
                    do x=1, nresi
                        if ((bsheet_assign(nresii + x) == 2) .or. (bsheet_assign(nresii + x) == 5)) anti = anti + 1
                    enddo                    
                else
                    do x=1, nresi
                        if ((bsheet_assign(nresii + x) == 2)) then
                            anti = anti + 1
                        elseif ((bsheet_assign(nresii + x) == 3)) then
                            para = para + 1
                        endif
                    enddo                    
                    
                    do x=1, nresj
                        if ((bsheet_assign(nresjj + x) == 2)) then
                            anti = anti + 1
                        elseif ((bsheet_assign(nresjj + x) == 3)) then
                            para = para + 1
                        endif
                    enddo
                    
                    !write(*,*) "error in determining direction in pair:", i, j, a, b
                    !write(*,*) "this can happen when a peptide has a turn"
                    !write(*,*) nh, co, nh2, co2
                    undetermined = undetermined + 1     
                endif 

                if (i == j .and. a <= 0) then
                    do m = 1, size(resNHi)
                        NH_i_pos = resNHi(m)
                        CO_j_pos = resCOj(m)
                        diff = abs(NH_i_pos - CO_j_pos)
                        found = .false.
                        CO_i_pos = -1
                        do k = 1, size(resNHi)
                            if (resNHi(k) == CO_j_pos) then
                                found = .true.
                                CO_i_pos = resCOj(k)
                                exit
                            endif
                        enddo
                        
                        if (found .and. CO_i_pos == NH_i_pos) then
                            if (diff >= 3 .and. diff <= 5) then
                                if ((NH_i_pos-2 >= 1) .and. (CO_j_pos+2 <= nresi)) then
                                    nh_basg = bsheet_assign(nresii + NH_i_pos-2)
                                    co_basg = bsheet_assign(nresii + CO_j_pos+2)
                                    if (((nh_basg == 2) .or. (nh_basg == 5)) .and. ((co_basg == 2) .or. (co_basg == 5))) then
                                        turn = turn + 1
                                    endif
                                endif
                            endif
                            
                        else if (.not. found) then 
                            if (diff >= 3 .and. diff <= 5) then
                                if ((NH_i_pos-2 >= 1) .and. (CO_j_pos+2 <= nresi)) then
                                    nh_basg = bsheet_assign(nresii + NH_i_pos-2)
                                    co_basg = bsheet_assign(nresii + CO_j_pos+2)
                                    if (((nh_basg == 2) .or. (nh_basg == 5)) .and. ((co_basg == 2) .or. (co_basg == 5))) then
                                        turn = turn + 2
                                    endif
                                endif
                            endif 
                        endif
                    enddo
                    
                endif
                
                deallocate(nh_i)
                deallocate(nh_j)
                deallocate(co_i)
                deallocate(co_j)
                
            endif
        enddo
    enddo

    counter = 0
    do k=1, total_residues
        if (bsheet_assign(k) >= 2) counter = counter + 1
    enddo
    content = real(counter)/real(total_residues)
    
    turn = turn / 2
    deallocate(bsheet_assign)
    end subroutine bsheet_content
    
    end module my_functions
    
    program bsheet_content_main
    
    use my_functions
    
    implicit none

    character(len=256)                      :: pep1, pep2, pep(10)
    character(len=100)                      :: energy_file, bptnr_file, pdb_file
    integer                                 :: anneal_steps, start_step, end_step
    integer                                 :: nc1, nc2, nb1, nb2, i, j, k, npts, nframes, n1, n2, nres1, nres2
    integer                                 :: nc(10), nb(10), nres(10), numseq
    integer                                 :: anti, para, undetermined, turn
    integer, allocatable                    :: bptnrs(:,:,:), hbmats(:,:,:), hb(:), anti_ratio(:), para_ratio(:)
    
    integer(8), allocatable                 :: collision_nums(:), contents(:), coll_series(:)
    integer(8)                              :: colltotal
    real(8)                                 :: ttotal, content
    real(8), allocatable                    :: t_series(:), t(:), temp(:)
    logical                                 :: flag1, flag2
    
    anneal_steps = -1
    start_step = -1
    end_step = -1

    do i = 0, 9999
        write(energy_file,'(A,I4.4,A)') 'results/run', i, '.energy'
        write(bptnr_file,'(A,I4.4,A)')  'results/run', i, '.bptnr'

        inquire(file=trim(energy_file), exist=flag1)
        inquire(file=trim(bptnr_file),  exist=flag2)

        if (flag1 .and. flag2) then
            if (start_step == -1) start_step = i
            end_step = i
        endif
    enddo

    if (start_step == -1) then
        write(*,'(A)') 'No matching run*.energy and run*.bptnr files were found.'
        stop
    endif

    write(*,'(A,I0)') 'First run: ', start_step
    write(*,'(A,I0)') 'Last run:  ', end_step
    
    
    ! Read peptide sequences and chain counts from the first available results/run*.pdb file.
    call analyze_pdb(pep, nc, nb, nres, numseq, pdb_file)

    if (numseq > 2) then
        write(*,'(A)') 'The PDB reader supports 10 peptide species, but the current beta-sheet analysis supports only two.'
        stop
    endif

    pep1 = pep(1)
    nc1 = nc(1)
    nb1 = nb(1)
    nres1 = nres(1)

    pep2 = ''
    nc2 = 0
    nb2 = 0
    nres2 = 0
    if (numseq == 2) then
        pep2 = pep(2)
        nc2 = nc(2)
        nb2 = nb(2)
        nres2 = nres(2)
    endif
    
    open(10,file='bsheet_content.txt', status='replace', action='write')
    close(10)

    !!!!!!!!!!!!!!!!! get start time !!!!!!!!!!!!!!!!!!!!!!!!!!
    ttotal = 0      ! Cumulative time before reading the current energy file
    colltotal = 0   ! Cumulative collisions before reading the current energy file
    if (start_step /= 0 .and. start_step /= 1) then
        do i = 0, start_step - 1
            write(energy_file, '(A,I4.4,A)') 'results/run', i, '.energy'
            call read_energy(energy_file, ttotal, colltotal, i, end_step, &
                anneal_steps, t_series, coll_series, npts)

            if (npts > 0) then
                ttotal =  t_series(npts)
            end if

            deallocate(t_series)
            deallocate(coll_series)
        end do
    end if
    
    k = 0
    do j = start_step, end_step
        write(bptnr_file, '(A,I4.4,A)') 'results/run', j, '.bptnr'
        write(energy_file, '(A,I4.4,A)') 'results/run', j, '.energy'
        
        call read_energy(energy_file, ttotal, colltotal, j, end_step, anneal_steps, &
            t_series, coll_series, npts)
        
        if (.not. allocated(t)) then
            allocate(t(size(t_series)))
            t = t_series
        else ! combine t_series (time in one simulaiton) to the whole time period
            n1 = size(t)
            n2 = size(t_series)
    
            allocate(temp(n1 + n2))
            if (n1 > 0) temp(1:n1) = t(1:n1)
            if (n2 > 0) temp(n1+1:n1+n2) = t_series(1:n2)
    
            deallocate(t)
            
            t = temp
            deallocate(temp)
        endif
        
        
        
        call read_bptnr_multi(bptnr_file, nc1, nc2, nb1, nb2, nframes, collision_nums, bptnrs, hbmats)

        if (nframes /= npts) then
            write(*,'(A,I0,A,I0,A,I0,A)') 'Run ', j, ' has ', npts, &
                ' retained energy rows but ', nframes, ' beta-sheet frames.'
            stop
        endif
    

        open(10,file='bsheet_content.txt', access="append")
        do i = 1, nframes
            
            call bsheet_content(bptnrs(:,:,i), hbmats(:,:,i), nc1, nc2, nb1, nb2, nres1, nres2, anti, para, turn, undetermined, content)
            write(10, "(F10.5,2X,F10.5,3I10,F15.5)") dble(coll_series(i))/1000000000.0, t_series(i), anti, para, turn, content

        enddo
        close(10)
        
        deallocate(t_series)
        deallocate(coll_series)
        deallocate(collision_nums)
        deallocate(bptnrs)
        deallocate(hbmats)
    enddo
    
    end program bsheet_content_main

