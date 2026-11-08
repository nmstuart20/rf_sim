//! Render blocks: half-open ranges of sample indices on one channel.

use std::iter::FusedIterator;

use crate::{SampleIndex, TimeError};

/// A half-open range `[start, end)` of sample indices on one channel.
///
/// Because `end` is exclusive and a u64, index `u64::MAX` is never inside a block.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct Block {
    start: SampleIndex,
    end: SampleIndex,
}

impl Block {
    /// The range `[start, end)`. [`TimeError::InvalidRange`] if `end < start`.
    pub fn new(start: SampleIndex, end: SampleIndex) -> Result<Self, TimeError> {
        if end < start {
            return Err(TimeError::InvalidRange);
        }
        Ok(Self { start, end })
    }

    /// First index in the block.
    pub const fn start(self) -> SampleIndex {
        self.start
    }

    /// One past the last index in the block.
    pub const fn end(self) -> SampleIndex {
        self.end
    }

    /// Number of samples.
    pub const fn len(self) -> u64 {
        self.end.get() - self.start.get()
    }

    /// Whether the block has no samples.
    pub const fn is_empty(self) -> bool {
        self.start.get() == self.end.get()
    }

    /// Whether `n` is inside the block.
    pub fn contains(self, n: SampleIndex) -> bool {
        self.start <= n && n < self.end
    }
}

/// The block grid of a channel: block `b` is `[b · len, (b + 1) · len)`.
///
/// Blocks are aligned to multiples of `len` from index 0, so a block's number does not
/// depend on where a run starts. The maximum length is pinned with the block
/// polynomials (prompt 2); here it only has to be nonzero.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub struct BlockSpec {
    len: u64,
}

impl BlockSpec {
    /// Blocks of `len` samples. [`TimeError::ZeroBlockLength`] if `len == 0`.
    pub fn new(len: u64) -> Result<Self, TimeError> {
        if len == 0 {
            return Err(TimeError::ZeroBlockLength);
        }
        Ok(Self { len })
    }

    /// Block length in samples.
    pub const fn block_len(self) -> u64 {
        self.len
    }

    /// Number of the block that contains `n`.
    pub const fn block_number_of(self, n: SampleIndex) -> u64 {
        n.get() / self.len
    }

    /// Block number `b`. [`TimeError::IndexOverflow`] if its end exceeds u64.
    pub fn block(self, b: u64) -> Result<Block, TimeError> {
        let start = b.checked_mul(self.len).ok_or(TimeError::IndexOverflow)?;
        let end = start
            .checked_add(self.len)
            .ok_or(TimeError::IndexOverflow)?;
        Ok(Block {
            start: SampleIndex::new(start),
            end: SampleIndex::new(end),
        })
    }

    /// Iterate over the grid blocks that intersect `span`, each clipped to `span`.
    ///
    /// The pieces are non-empty, contiguous and in order, and together cover `span`
    /// exactly. Only the first and last piece can be shorter than `len`.
    pub const fn tile(self, span: Block) -> Blocks {
        Blocks {
            spec: self,
            next: span.start.get(),
            end: span.end.get(),
        }
    }

    /// Iterate over whole grid blocks starting at block `first`, until the u64 index
    /// space ends (the final piece is clipped at `u64::MAX`).
    pub fn iter_from(self, first: u64) -> Result<Blocks, TimeError> {
        let start = first
            .checked_mul(self.len)
            .ok_or(TimeError::IndexOverflow)?;
        Ok(Blocks {
            spec: self,
            next: start,
            end: u64::MAX,
        })
    }
}

/// Iterator over block pieces, created by [`BlockSpec::tile`] or [`BlockSpec::iter_from`].
/// Yields `(block number, piece)`.
#[derive(Clone, Debug)]
pub struct Blocks {
    spec: BlockSpec,
    next: u64,
    end: u64,
}

impl Iterator for Blocks {
    type Item = (u64, Block);

    fn next(&mut self) -> Option<Self::Item> {
        if self.next >= self.end {
            return None;
        }
        let number = self.next / self.spec.len;
        // The grid end can exceed u64 only for the last block; clip it to `end`.
        let grid_end = number
            .checked_add(1)
            .and_then(|b| b.checked_mul(self.spec.len))
            .unwrap_or(u64::MAX);
        let piece_end = grid_end.min(self.end);
        let piece = Block {
            start: SampleIndex::new(self.next),
            end: SampleIndex::new(piece_end),
        };
        self.next = piece_end;
        Some((number, piece))
    }
}

impl FusedIterator for Blocks {}

#[cfg(test)]
mod tests {
    use super::*;

    fn span(a: u64, b: u64) -> Block {
        Block::new(SampleIndex::new(a), SampleIndex::new(b)).unwrap()
    }

    #[test]
    fn tiles_with_partial_ends() {
        let spec = BlockSpec::new(10).unwrap();
        let pieces: Vec<_> = spec
            .tile(span(5, 32))
            .map(|(b, p)| (b, p.start().get(), p.end().get()))
            .collect();
        assert_eq!(pieces, [(0, 5, 10), (1, 10, 20), (2, 20, 30), (3, 30, 32)]);
        assert_eq!(spec.tile(span(7, 7)).count(), 0);
    }

    #[test]
    fn rejects_bad_input() {
        assert_eq!(BlockSpec::new(0), Err(TimeError::ZeroBlockLength));
        assert_eq!(
            Block::new(SampleIndex::new(2), SampleIndex::new(1)),
            Err(TimeError::InvalidRange)
        );
        let spec = BlockSpec::new(1 << 40).unwrap();
        assert_eq!(spec.block(1 << 24), Err(TimeError::IndexOverflow));
        assert!(spec.iter_from(1 << 24).is_err());
    }

    #[test]
    fn iter_from_reaches_end_of_index_space() {
        let spec = BlockSpec::new(1 << 62).unwrap();
        let pieces: Vec<_> = spec.iter_from(0).unwrap().collect();
        assert_eq!(pieces.len(), 4);
        assert_eq!(pieces[3].1.end().get(), u64::MAX);
        assert_eq!(spec.block(2).unwrap(), pieces[2].1);
    }

    #[test]
    fn one_ms_blocks_at_61_44_msps() {
        let spec = BlockSpec::new(61_440).unwrap();
        let b = spec.block(3).unwrap();
        assert_eq!(
            (b.start().get(), b.end().get(), b.len()),
            (184_320, 245_760, 61_440)
        );
        assert!(b.contains(SampleIndex::new(184_320)));
        assert!(!b.contains(SampleIndex::new(245_760)));
        assert_eq!(spec.block_number_of(SampleIndex::new(245_759)), 3);
    }
}
