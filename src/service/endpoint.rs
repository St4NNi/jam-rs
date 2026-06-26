//! The QUIC endpoint of a server or client inside one cluster: no relays and no address lookup;
//! peers are dialled by node id and address.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::ServiceError;
use iroh::endpoint::presets;
use iroh::{Endpoint, EndpointAddr, EndpointId, RelayMode, SecretKey};
use std::net::SocketAddr;

/// Binds an endpoint whose node id is the public key of `key`. Without `addresses` it binds
/// free ports on all interfaces.
pub async fn bind(key: SecretKey, addresses: &[SocketAddr]) -> Result<Endpoint, ServiceError> {
    let mut builder = Endpoint::builder(presets::Minimal)
        .relay_mode(RelayMode::Disabled)
        .clear_address_lookup()
        .secret_key(key);
    for address in addresses {
        builder = builder.bind_addr(*address)?;
    }
    Ok(builder.bind().await?)
}

/// Reads `NODE_ID@ADDRESS`, the form `jam serve` prints, such as `4f2a...@10.0.0.7:4100`.
pub fn parse_address(text: &str) -> Result<EndpointAddr, ServiceError> {
    let bad = || ServiceError::Address(text.to_string());
    let (id, address) = text.split_once('@').ok_or_else(bad)?;
    let id: EndpointId = id.parse().map_err(|_| bad())?;
    let address: SocketAddr = address.parse().map_err(|_| bad())?;
    Ok(EndpointAddr::new(id).with_ip_addr(address))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn address_reads_back() {
        let id = SecretKey::from_bytes(&[7; 32]).public();
        let text = format!("{id}@127.0.0.1:4100");
        let address = parse_address(&text).unwrap();
        assert_eq!(address.id, id);
        let ip: Vec<_> = address.ip_addrs().collect();
        assert_eq!(ip, [&"127.0.0.1:4100".parse::<SocketAddr>().unwrap()]);
        for bad in [
            "127.0.0.1:4100",
            "abc@127.0.0.1:4100",
            &format!("{id}@host"),
        ] {
            assert!(parse_address(bad).is_err(), "{bad}");
        }
    }
}
