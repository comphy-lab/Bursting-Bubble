/**
# adapt_wavelet_limited.h — adapt_wavelet with a position-dependent max level

Regional refinement ceiling: `MLFun(x,y,z)` returns the allowed maximum
level for each cell, replacing the single global `maxlevel` argument of
`adapt_wavelet`.

Provenance: the idea and the original implementation are C. Pairetti's
Basilisk sandbox utility
(http://basilisk.fr/sandbox/pairetti/bag_mode/adapt_wavelet_limited.h,
fetched 2026-07-02, byte-identical to the basilisk-C wiki mirror). That
file targets an older Basilisk API (`s.coarsen`/`no_coarsen`, pre-block
fields) and no longer compiles against the current tree. This version
re-applies Pairetti's modification — `cellMAX = MLFun(x,y,z)` evaluated
inside `foreach_cell()`, substituted for `maxlevel` in the three level
comparisons — onto the CURRENT `adapt_wavelet` from
`$BASILISK/grid/tree-common.h` (the copy this project compiles against,
2026-07-02), keeping everything else byte-for-byte.

Used by the drill solver to hold full resolution on the jet (at/above the
base) while capping the satellite / floor-remnant zone below the base —
repeated pinch singularities there must not be resolved at full depth
(case-1009 lesson).
*/

trace
astats adapt_wavelet_limited (scalar * slist,     // list of scalars
		      double * max,               // tolerance for each scalar
		      int (*MLFun)(double,double,double), // max level from position
		      int minlevel = 1,           // minimum level of refinement
		      scalar * list = all)        // list of fields to update
{
  scalar * ilist = list;

  if (is_constant(cm)) {
    if (list == NULL || list == all)
      list = list_copy (all);
    boundary (list);
    restriction (slist);
  }
  else {
    if (list == NULL || list == all) {
      list = list_copy ({cm, fm});
      for (scalar s in all)
	list = list_add (list, s);
    }
    boundary (list);
    scalar * listr = list_concat (slist, {cm});
    restriction (listr);
    free (listr);
  }

  astats st = {0, 0};
  scalar * listc = NULL;
  for (scalar s in list)
    listc = list_add_depend (listc, s);

  // refinement
  if (minlevel < 1)
    minlevel = 1;
  tree->refined.n = 0;
  static const int refined = 1 << user, too_fine = 1 << (user + 1);
  foreach_cell() {
    int cellMAX = MLFun (x, y, z);   // position-dependent ceiling (Pairetti)
    if (is_active(cell)) {
      static const int too_coarse = 1 << (user + 2);
      if (is_leaf (cell)) {
	if (cell.flags & too_coarse) {
	  cell.flags &= ~too_coarse;
	  refine_cell (point, listc, refined, &tree->refined);
	  st.nf++;
	}
	continue;
      }
      else { // !is_leaf (cell)
	if (cell.flags & refined) {
	  // cell has already been refined, skip its children
	  cell.flags &= ~too_coarse;
	  continue;
	}
	// check whether the cell or any of its children is local
	bool local = is_local(cell);
	if (!local)
	  foreach_child()
	    if (is_local(cell)) {
	      local = true; break;
	    }
	if (local) {
	  int i = 0;
	  static const int just_fine = 1 << (user + 3);
	  for (scalar s in slist) {
	    double emax = max[i++], sc[(1 << dimension)*s.block];
	    double * b = sc;
	    foreach_child()
	      foreach_blockf(s)
	        *b++ = s[];
	    s.prolongation (point, s);
	    b = sc;
	    foreach_child()
	      foreach_blockf(s) {
	        double e = fabs(*b - s[]);
		if (e > emax && level < cellMAX) {
		  cell.flags &= ~too_fine;
		  cell.flags |= too_coarse;
		}
		else if ((e <= emax/1.5 || level > cellMAX) &&
			 !(cell.flags & (too_coarse|just_fine))) {
		  if (level >= minlevel)
		    cell.flags |= too_fine;
		}
		else if (!(cell.flags & too_coarse)) {
		  cell.flags &= ~too_fine;
		  cell.flags |= just_fine;
		}
		s[] = *b++;
	      }
	  }
	  foreach_child() {
	    cell.flags &= ~just_fine;
	    if (!is_leaf(cell)) {
	      cell.flags &= ~too_coarse;
	      if (level >= cellMAX)
		cell.flags |= too_fine;
	    }
	    else if (!is_active(cell))
	      cell.flags &= ~too_coarse;
	  }
	}
      }
    }
    else // inactive cell
      continue;
  }
  mpi_boundary_refine (listc);

  // coarsening
  // the loop below is only necessary to ensure symmetry of 2:1 constraint
  for (int l = depth(); l >= 0; l--) {
    foreach_cell()
      if (!is_boundary(cell)) {
	if (level == l) {
	  if (!is_leaf(cell)) {
	    if (cell.flags & refined)
	      // cell was refined previously, unset the flag
	      cell.flags &= ~(refined|too_fine);
	    else if (cell.flags & too_fine) {
	      if (is_local(cell) && coarsen_cell (point, listc))
		st.nc++;
	      cell.flags &= ~too_fine; // do not coarsen parent
	    }
	  }
	  if (cell.flags & too_fine)
	    cell.flags &= ~too_fine;
	  else if (level > 0 && (aparent(0).flags & too_fine))
	    aparent(0).flags &= ~too_fine;
	  continue;
	}
	else if (is_leaf(cell))
	  continue;
      }
    mpi_boundary_coarsen (l, too_fine);
  }
  free (listc);

  mpi_all_reduce (st.nf, MPI_INT, MPI_SUM);
  mpi_all_reduce (st.nc, MPI_INT, MPI_SUM);
  if (st.nc || st.nf)
    mpi_boundary_update (list);

  if (list != ilist)
    free (list);

  return st;
}
