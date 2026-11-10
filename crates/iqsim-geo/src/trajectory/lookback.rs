//! Bounded backward Taylor extrapolation for retarded trajectory evaluation.

use std::sync::Arc;

use iqsim_core::{AcqInstant, TruthDuration, TruthTime};

use super::{Span, State, Trajectory, seconds_since};
use crate::{Ecef, GeoError};

/// Largest permitted trajectory lookback: one second.
pub const MAX_LOOKBACK: TruthDuration = TruthDuration::from_ps(1_000_000_000_000);

/// A trajectory with a bounded, second-order extrapolation before its start.
///
/// Translation is C2 at the join. Attitude and body angular rate are held at their
/// initial values throughout extrapolation. An unbounded-start trajectory passes
/// through unchanged, including its span and knots.
#[derive(Clone, Debug)]
pub struct Lookback {
    inner: Arc<dyn Trajectory>,
    span: Span,
    initial: Option<(TruthTime, State)>,
}

impl Lookback {
    /// Wrap `inner` with a lookback in `(0, MAX_LOOKBACK]`.
    ///
    /// Rejects invalid durations and overflow of the extended span start. For a
    /// bounded start, evaluates and caches its state; any trajectory error propagates.
    pub fn new(inner: Arc<dyn Trajectory>, lookback: TruthDuration) -> Result<Self, GeoError> {
        if lookback <= TruthDuration::ZERO || lookback > MAX_LOOKBACK {
            return Err(GeoError::Trajectory("lookback must be in (0, 1 s]".into()));
        }
        let mut span = inner.span();
        let initial = if let Some(start) = span.start {
            span.start = Some(start.checked_sub(lookback).map_err(|_| {
                GeoError::Trajectory("lookback span start overflows truth time".into())
            })?);
            Some((start, inner.state_at(start)?))
        } else {
            None
        };
        Ok(Self {
            inner,
            span,
            initial,
        })
    }
}

impl Trajectory for Lookback {
    fn state(&self, t: AcqInstant) -> Result<State, GeoError> {
        let Some((start, initial)) = self.initial else {
            return self.inner.state(t);
        };
        if !self.span.contains(t) {
            return Err(GeoError::OutsideSpan {
                t_ps: t.ps().as_ps(),
                start_ps: self.span.start.unwrap_or(TruthTime::MIN).as_ps(),
                end_ps: self.span.end.unwrap_or(TruthTime::MAX).as_ps(),
            });
        }
        if t >= AcqInstant::from_truth_time(start) {
            return self.inner.state(t);
        }
        let dt = seconds_since(t, start);
        // Fixed component-wise order: p + (v*dt + (a*(0.5*dt))*dt).
        let half_dt = 0.5 * dt;
        let quadratic = (initial.acceleration * half_dt) * dt;
        let position = initial.position.vec() + (initial.velocity * dt + quadratic);
        Ok(State {
            position: Ecef::from_vec(position),
            velocity: initial.velocity + initial.acceleration * dt,
            ..initial
        })
    }

    fn span(&self) -> Span {
        self.span
    }

    fn knots(&self, from: TruthTime, to: TruthTime) -> Vec<TruthTime> {
        let mut knots = self.inner.knots(from, to);
        if let Some((start, _)) = self.initial {
            if from <= start && start <= to {
                knots.push(start);
            }
            knots.sort_unstable();
            knots.dedup();
        }
        knots
    }
}
