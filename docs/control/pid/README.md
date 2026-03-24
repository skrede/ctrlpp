# PID

A policy-based PID controller with compile-time feature composition. The template
signature `pid<Scalar, NX, NU, NY, Policies...>` lets you opt into exactly the
features you need -- anti-windup, derivative filtering, setpoint weighting,
velocity form, ISA form, feed-forward, rate limiting, and performance assessment --
without paying for what you don't use.

The base PID provides proportional, integral, and derivative action with output
clamping. Each policy extends the controller with one additional behaviour,
composed at compile time via variadic template packs.

## Pages

- [anti-windup.md](anti-windup.md) -- Back-calculation, clamping, and conditional integration
- [derivative-filter.md](derivative-filter.md) -- Low-pass filtered derivative term
- [setpoint-filter.md](setpoint-filter.md) -- Setpoint weighting and filtering
- [velocity-form.md](velocity-form.md) -- Incremental (velocity) form PID
- [isa-form.md](isa-form.md) -- ISA standard form (Ti, Td parametrisation)
- [feed-forward.md](feed-forward.md) -- Additive feed-forward policy
- [rate-limit.md](rate-limit.md) -- Control signal rate limiting
- [performance.md](performance.md) -- IAE, ISE, ITAE, and oscillation detection

## Composition Guide

See [PID Composition](../../guides/pid/composition.md) for a walkthrough of
combining multiple policies into a single controller.

## Theory

For the mathematical background, see [PID Theory](../../reference/pid-theory.md).
