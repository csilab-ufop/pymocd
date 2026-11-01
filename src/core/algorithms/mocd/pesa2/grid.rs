//! The PESA-II hyper-grid: niche binning, squeeze factors, region selection.
//! This Source Code Form is subject to the terms of The GNU General Public License v3.0
//! Copyright 2025 - Guilherme Santos. If a copy of the MPL was not distributed with this
//! file, You can obtain one at https://www.gnu.org/licenses/gpl-3.0.html

use rand::RngExt;

use super::solution::Member;

/// Hyper-grid resolution per objective axis: `NICHSIZE` of the authors'
/// implementation, a 10 x 10 grid for the two objectives.
pub const GRID_DIVISIONS: usize = 10;

/// Each axis is normalised to the range the members *currently* span, so the
/// grid moves with the archive and must be rebuilt after any change to it.
/// `occ`, the per-cell occupancy, is the squeeze factor vector.
pub fn assign_cells(members: &mut [Member], divisions: usize, occ: &mut Vec<usize>) {
    let obj_len = members[0].solution.objectives.len();
    let mut min_v = vec![f64::INFINITY; obj_len];
    let mut max_v = vec![f64::NEG_INFINITY; obj_len];
    for m in members.iter() {
        for k in 0..obj_len {
            let v = m.solution.objectives[k];
            min_v[k] = min_v[k].min(v);
            max_v[k] = max_v[k].max(v);
        }
    }

    for m in members.iter_mut() {
        let mut cell = 0usize;
        for k in 0..obj_len {
            let span = max_v[k] - min_v[k];
            let norm = if span.abs() < f64::EPSILON {
                0.0
            } else {
                (m.solution.objectives[k] - min_v[k]) / span
            };
            let bin = ((norm * divisions as f64) as usize).min(divisions - 1);
            cell = cell * divisions + bin;
        }
        m.cell = cell;
    }

    occ.clear();
    occ.resize(divisions.pow(obj_len as u32), 0);
    for m in members.iter() {
        occ[m.cell] += 1;
    }
}

/// Region-based selection as the authors' implementation performs it: a
/// uniformly random *occupied* cell, then a uniformly random member of it. A
/// sparse cell is as likely to be drawn as a crowded one, which spreads the
/// parents along the front.
pub fn select_member<'a>(
    ep: &'a [Member],
    occ: &[usize],
    occupied: &[usize],
    rng: &mut impl rand::Rng,
) -> &'a Member {
    let cell = occupied[rng.random_range(0..occupied.len())];
    let k = rng.random_range(0..occ[cell]);
    ep.iter()
        .filter(|m| m.cell == cell)
        .nth(k)
        .expect("occupancy out of sync with EP")
}
