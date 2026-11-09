//! Lint diagnostics with source locations.

use std::fmt::Write as _;
use std::ops::Range;

/// One problem found in a scenario, located by byte range in its source.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Diagnostic {
    /// Byte range of the offending text, if known.
    pub span: Option<Range<usize>>,
    /// What is wrong.
    pub message: String,
    /// A related location, such as where a duplicate name was first defined.
    pub note: Option<(Range<usize>, String)>,
}

impl Diagnostic {
    /// A diagnostic at `span`.
    pub fn new(span: Range<usize>, message: impl Into<String>) -> Self {
        Self {
            span: Some(span),
            message: message.into(),
            note: None,
        }
    }

    /// Adds a related location.
    pub fn with_note(mut self, span: Range<usize>, note: impl Into<String>) -> Self {
        self.note = Some((span, note.into()));
        self
    }

    /// The diagnostic for a TOML parse or type error.
    pub fn from_toml(error: &toml::de::Error) -> Self {
        Self {
            span: error.span(),
            message: error.message().trim_end().to_owned(),
            note: None,
        }
    }

    /// 1-based line and column of the start of the span, if any.
    pub fn location(&self, source: &str) -> Option<(usize, usize)> {
        self.span.as_ref().map(|s| line_col(source, s.start))
    }

    /// Renders `path:line:col: error: message`, the source line and a caret under
    /// the span, then the note in the same form.
    pub fn render(&self, path: &str, source: &str) -> String {
        let mut out = String::new();
        render_one(
            &mut out,
            path,
            source,
            self.span.as_ref(),
            "error",
            &self.message,
        );
        if let Some((span, note)) = &self.note {
            render_one(&mut out, path, source, Some(span), "note", note);
        }
        out
    }
}

fn render_one(
    out: &mut String,
    path: &str,
    source: &str,
    span: Option<&Range<usize>>,
    level: &str,
    message: &str,
) {
    let Some(span) = span else {
        let _ = writeln!(out, "{path}: {level}: {message}");
        return;
    };
    let (line, col) = line_col(source, span.start);
    let _ = writeln!(out, "{path}:{line}:{col}: {level}: {message}");
    let start = clamp_to_char(source, span.start);
    let line_start = source[..start].rfind('\n').map_or(0, |i| i + 1);
    let line_end = source[start..]
        .find('\n')
        .map_or(source.len(), |i| start + i);
    let text = source[line_start..line_end].trim_end_matches('\r');
    let end = clamp_to_char(source, span.end.clamp(start, line_end));
    let width = source[start..end].chars().count().max(1);
    let gutter = line.to_string().len();
    let _ = writeln!(out, "{line:>gutter$} | {text}");
    let _ = writeln!(
        out,
        "{:>gutter$} | {}{}",
        "",
        " ".repeat(col - 1),
        "^".repeat(width)
    );
}

/// The largest char boundary at or before `offset` (clamped to the source).
fn clamp_to_char(source: &str, offset: usize) -> usize {
    let mut i = offset.min(source.len());
    while !source.is_char_boundary(i) {
        i -= 1;
    }
    i
}

/// 1-based line and column (in characters) of byte `offset` in `source`.
pub fn line_col(source: &str, offset: usize) -> (usize, usize) {
    let offset = clamp_to_char(source, offset);
    let before = &source[..offset];
    let line = before.matches('\n').count() + 1;
    let line_start = before.rfind('\n').map_or(0, |i| i + 1);
    (line, before[line_start..].chars().count() + 1)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn locations_are_one_based_characters() {
        let src = "a = 1\nµ = \"x\"\n";
        assert_eq!(line_col(src, 0), (1, 1));
        assert_eq!(line_col(src, 4), (1, 5));
        assert_eq!(line_col(src, 6), (2, 1));
        // 'µ' is two bytes, one column.
        assert_eq!(line_col(src, 8), (2, 2));
        assert_eq!(line_col(src, 999), (3, 1));
    }

    #[test]
    fn renders_caret_under_span() {
        let src = "x = 1\nfreq = \"1090\"\n";
        let d = Diagnostic::new(13..19, "frequency has no unit").with_note(0..1, "see here");
        assert_eq!(
            d.render("s.toml", src),
            "s.toml:2:8: error: frequency has no unit\n\
             2 | freq = \"1090\"\n  |        ^^^^^^\n\
             s.toml:1:1: note: see here\n\
             1 | x = 1\n  | ^\n"
        );
        let d = Diagnostic {
            span: None,
            message: "m".into(),
            note: None,
        };
        assert_eq!(d.render("s.toml", src), "s.toml: error: m\n");
    }
}
