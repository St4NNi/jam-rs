use anyhow::{Context, Result, ensure};
use jam_rs::shared_reader::SharedReader;
use jam_rs::shared_seed::SharedKey;
use std::path::PathBuf;

fn main() -> Result<()> {
    let args = std::env::args_os()
        .skip(1)
        .map(PathBuf::from)
        .collect::<Vec<_>>();
    ensure!(
        args.len() == 3,
        "usage: shared_compare BASELINE PACKED KEYS_JSON"
    );
    let baseline = SharedReader::open(&args[0])?;
    let packed = SharedReader::open(&args[1])?;
    baseline.verify_checksum()?;
    packed.verify_checksum()?;
    ensure!(baseline.core_count() == packed.core_count());
    ensure!(baseline.occurrence_count() == packed.occurrence_count());
    ensure!(baseline.manifest_sha256() == packed.manifest_sha256());
    let keys: Vec<(u32, u32, u8)> = serde_json::from_reader(std::fs::File::open(&args[2])?)?;
    ensure!(keys.len() <= 65_536, "representative key budget");
    let mut groups = 0u64;
    let mut memberships = 0u64;
    let mut placements = 0u64;
    for &(core, context, length) in &keys {
        let key = SharedKey {
            core,
            context,
            length,
        };
        let left = baseline.find(key)?;
        let right = packed.find(key)?;
        ensure!(left.is_some() == right.is_some(), "group presence");
        let Some(left) = left else { continue };
        let right = right.context("packed group")?;
        ensure!(
            left.member_count() == right.member_count(),
            "group member count"
        );
        ensure!(
            left.occurrence_count() == right.occurrence_count(),
            "group occurrence count"
        );
        let members = baseline.members(left)?;
        let packed_members = packed.members(right)?;
        ensure!(members.len() == packed_members.len(), "memberships");
        for (member, other) in members.into_iter().zip(packed_members) {
            ensure!(
                member.metagenome_id == other.metagenome_id,
                "member identity"
            );
            ensure!(
                member.occurrence_count() == other.occurrence_count(),
                "member count"
            );
            ensure!(
                packed.member(right, member.metagenome_id)? == Some(other),
                "direct member"
            );
            let mut at = 0;
            while at < member.occurrence_count() {
                let expected = baseline.occurrence_block(left, member, at, 4096)?;
                let actual = packed.occurrence_block(right, other, at, 4096)?;
                ensure!(actual == expected && !actual.is_empty(), "placement block");
                at += actual.len() as u64;
                placements += actual.len() as u64;
            }
            memberships += 1;
        }
        groups += 1;
    }
    println!(
        "{}",
        serde_json::json!({
            "status": "equal", "requested_keys": keys.len(), "groups": groups,
            "memberships": memberships, "context_associated_placements": placements,
            "scope": "deterministic representative keys; full population streaming audit separate"
        })
    );
    Ok(())
}
