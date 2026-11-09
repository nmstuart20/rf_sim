//! Names and references, with their source spans.

use std::borrow::Cow;
use std::ops::Range;

use schemars::JsonSchema;
use serde::{Deserialize, Deserializer};
use toml::Spanned;

/// Checks the identifier charset; returns the reason it is invalid.
fn check(text: &str, allow_path: bool) -> Result<(), String> {
    let mut chars = text.chars();
    let Some(first) = chars.next() else {
        return Err("a name cannot be empty".to_owned());
    };
    if !first.is_ascii_alphanumeric() {
        return Err(format!(
            "\"{text}\" must start with an ASCII letter or digit"
        ));
    }
    let extra: &[char] = if allow_path {
        &['_', '-', '.', '/']
    } else {
        &['_', '-']
    };
    if let Some(bad) = text
        .chars()
        .find(|c| !c.is_ascii_alphanumeric() && !extra.contains(c))
    {
        let allowed = if allow_path {
            "letters, digits, '_', '-', '.' and '/'"
        } else {
            "letters, digits, '_' and '-'"
        };
        return Err(format!(
            "\"{text}\" contains '{bad}'; names use ASCII {allowed}"
        ));
    }
    Ok(())
}

macro_rules! spanned_name {
    ($(#[$doc:meta])* $ty:ident, $allow_path:literal, $pattern:literal, $description:literal) => {
        $(#[$doc])*
        #[derive(Clone, Debug)]
        pub struct $ty(Spanned<String>);

        impl $ty {
            /// The name text.
            pub fn as_str(&self) -> &str {
                self.0.get_ref()
            }

            /// Byte range of the value in the scenario source.
            pub fn span(&self) -> Range<usize> {
                self.0.span()
            }
        }

        impl PartialEq for $ty {
            /// Names compare by text only, not by where they were written.
            fn eq(&self, other: &Self) -> bool {
                self.as_str() == other.as_str()
            }
        }

        impl Eq for $ty {}

        impl<'de> Deserialize<'de> for $ty {
            fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
                let spanned = Spanned::<String>::deserialize(d)?;
                check(spanned.get_ref(), $allow_path).map_err(serde::de::Error::custom)?;
                Ok(Self(spanned))
            }
        }

        impl JsonSchema for $ty {
            fn schema_name() -> Cow<'static, str> {
                stringify!($ty).into()
            }

            fn json_schema(_: &mut schemars::SchemaGenerator) -> schemars::Schema {
                schemars::json_schema!({
                    "type": "string",
                    "pattern": $pattern,
                    "description": $description,
                })
            }
        }
    };
}

spanned_name!(
    /// The name of a scenario entity (platform, emitter, receiver, channel, mount,
    /// antenna port, sink or resource), or a reference to one. ASCII letters,
    /// digits, `_` and `-`, starting with a letter or digit.
    Ident,
    false,
    "^[A-Za-z0-9][A-Za-z0-9_-]*$",
    "An entity name or a reference to one: ASCII letters, digits, '_' and '-'."
);

spanned_name!(
    /// A reference to a library item (antenna, waveform, behavior), found by name on
    /// the library search path. Like [`Ident`], and may also contain `.` and `/`.
    LibraryRef,
    true,
    "^[A-Za-z0-9][A-Za-z0-9_./-]*$",
    "A library item name, resolved on the library search path: ASCII letters, digits, '_', '-', '.' and '/'."
);

#[cfg(test)]
mod tests {
    use super::check;

    #[test]
    fn charset() {
        assert!(check("rx-1_a", false).is_ok());
        assert!(check("9a", false).is_ok());
        assert!(check("", false).is_err());
        assert!(check("_a", false).is_err());
        assert!(check("a b", false).unwrap_err().contains("' '"));
        assert!(check("ssr/mode_ac.v1", false).is_err());
        assert!(check("ssr/mode_ac.v1", true).is_ok());
    }
}
