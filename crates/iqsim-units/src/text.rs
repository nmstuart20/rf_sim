//! Shared parsing, serde and JSON Schema support for quantities written as strings.

use std::fmt;
use std::marker::PhantomData;

use serde::de::{self, Visitor};

use crate::decimal::{Decimal, split_number};
use crate::error::{QuantityKind, UnitError};

/// ECMA-262 regex for the number part of a quantity (matches `decimal::number_len`).
pub(crate) const NUMBER_PATTERN: &str =
    r"[-+]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][-+]?[0-9]+)?";

/// A quantity parsed from a string with a unit.
pub(crate) trait TextQuantity: Sized {
    const KIND: &'static QuantityKind;
    fn parse_text(text: &str) -> Result<Self, UnitError>;
}

/// Splits `text` into an exact number and a non-empty unit.
pub(crate) fn number_and_unit<'a>(
    kind: &'static QuantityKind,
    text: &'a str,
) -> Result<(Decimal, &'a str), UnitError> {
    let trimmed = text.trim();
    if trimmed.is_empty() {
        return Err(UnitError::Empty { kind });
    }
    let (number, unit) = split_number(trimmed);
    let value = Decimal::parse(number).ok_or_else(|| UnitError::InvalidNumber {
        kind,
        text: text.to_owned(),
    })?;
    if unit.is_empty() {
        return Err(UnitError::MissingUnit {
            kind,
            text: text.to_owned(),
        });
    }
    Ok((value, unit))
}

/// Rejects non-finite float results (overflowing exponents).
pub(crate) fn finite(
    kind: &'static QuantityKind,
    text: &str,
    value: f64,
) -> Result<f64, UnitError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(UnitError::OutOfRange {
            kind,
            text: text.to_owned(),
            limit: "the value overflows a 64-bit float".to_owned(),
        })
    }
}

/// `^\s*NUMBER\s*(?:unit|unit...)\s*$`.
pub(crate) fn unit_pattern(units: &[&str]) -> String {
    format!(r"^\s*{NUMBER_PATTERN}\s*(?:{})\s*$", units.join("|"))
}

/// The schema description for a quantity kind.
pub(crate) fn description(kind: &QuantityKind) -> String {
    format!(
        "A {} with a unit ({}), written as a string, e.g. \"{}\". Units are case-sensitive.",
        kind.name,
        kind.unit_list(),
        kind.example
    )
}

/// Serde visitor that parses strings and explains bare numbers.
pub(crate) struct TextVisitor<T>(PhantomData<T>);

impl<T> TextVisitor<T> {
    pub(crate) const fn new() -> Self {
        Self(PhantomData)
    }
}

impl<T: TextQuantity> Visitor<'_> for TextVisitor<T> {
    type Value = T;

    fn expecting(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "a {} string such as \"{}\"",
            T::KIND.name,
            T::KIND.example
        )
    }

    fn visit_str<E: de::Error>(self, v: &str) -> Result<T, E> {
        T::parse_text(v).map_err(E::custom)
    }

    fn visit_i64<E: de::Error>(self, v: i64) -> Result<T, E> {
        Err(bare(T::KIND, v.to_string()))
    }

    fn visit_u64<E: de::Error>(self, v: u64) -> Result<T, E> {
        Err(bare(T::KIND, v.to_string()))
    }

    fn visit_f64<E: de::Error>(self, v: f64) -> Result<T, E> {
        Err(bare(T::KIND, v.to_string()))
    }
}

fn bare<E: de::Error>(kind: &'static QuantityKind, number: String) -> E {
    E::custom(UnitError::BareNumber { kind, number })
}

/// Implements `FromStr`, `Deserialize` and `JsonSchema` for a [`TextQuantity`].
/// `$pattern` is an expression giving the schema regex, or `None`.
macro_rules! text_quantity_impls {
    ($ty:ident, $pattern:expr) => {
        impl ::std::str::FromStr for $ty {
            type Err = $crate::UnitError;
            fn from_str(text: &str) -> Result<Self, Self::Err> {
                <Self as $crate::text::TextQuantity>::parse_text(text)
            }
        }

        impl<'de> ::serde::Deserialize<'de> for $ty {
            fn deserialize<D: ::serde::Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
                d.deserialize_any($crate::text::TextVisitor::<Self>::new())
            }
        }

        impl ::schemars::JsonSchema for $ty {
            fn schema_name() -> ::std::borrow::Cow<'static, str> {
                stringify!($ty).into()
            }

            fn json_schema(_: &mut ::schemars::SchemaGenerator) -> ::schemars::Schema {
                let kind = <Self as $crate::text::TextQuantity>::KIND;
                let pattern: Option<String> = $pattern;
                let mut schema = ::schemars::json_schema!({
                    "type": "string",
                    "description": $crate::text::description(kind),
                    "examples": [kind.example],
                });
                if let Some(pattern) = pattern {
                    schema.insert("pattern".to_owned(), pattern.into());
                }
                schema
            }
        }
    };
}

pub(crate) use text_quantity_impls;
