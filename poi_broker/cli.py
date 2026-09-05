"""Flask CLI commands."""

import click
from flask import Flask


def register_cli(app: Flask) -> None:
    """Register one-off operator commands on ``app``."""

    @app.cli.command('fail-stale-exports')
    @click.option(
        '--force',
        is_flag=True,
        help=(
            'Fail all PENDING/RUNNING exports immediately, ignoring age. '
            'Use after Ctrl+C in memory mode, or when the Huey consumer is down. '
            'Aborts any in-flight export.'
        ),
    )
    def fail_stale_exports(force: bool) -> None:
        """Fail stuck PENDING/RUNNING export tasks so a user can retry."""
        from .tasks import reset_stale_export_tasks

        n = reset_stale_export_tasks(ignore_age=force)
        click.echo(f'Marked {n} export task(s) FAILED')
