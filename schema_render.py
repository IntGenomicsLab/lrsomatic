"""
schema_render.py
MkDocs macros module for rendering nextflow_schema.json as HTML parameter docs.
Loaded by mkdocs-macros-plugin via `module_name: schema_render` in mkdocs.yml.
"""

import json
import os
import html as html_module

DEFAULT_ICON = "fas fa-cog"


def define_env(env):
    env.variables["git_ref"] = os.environ.get("DOCS_REF", "dev")

    @env.macro
    def render_schema(
        schema_path=None,
        show_hidden=False,
        section_title=None,
    ):
        """
        Render nextflow_schema.json as styled HTML parameter tables.

        Args:
            schema_path:   Path to schema JSON. Defaults to nextflow_schema.json at
                           the project root (directory containing mkdocs.yml).
            show_hidden:   If True, renders only parameters with hidden=true.
                           If False (default), renders only visible parameters.
            section_title: Optional H2 heading rendered above all groups.

        Usage in markdown:
            {{ render_schema() }}
            {{ render_schema(show_hidden=True, section_title="Hidden / Advanced Parameters") }}
        """

        if schema_path is None:
            schema_path = os.path.join(env.project_dir, "nextflow_schema.json")

        with open(schema_path, "r", encoding="utf-8") as f:
            schema = json.load(f)

        # Support both $defs (JSON Schema 2020-12) and legacy definitions
        defs = schema.get("$defs", schema.get("definitions", {}))

        # Use allOf to respect the canonical group ordering defined in the schema
        all_of = schema.get("allOf", [])
        if all_of:
            ordered_group_keys = [
                ref["$ref"].split("/")[-1]
                for ref in all_of
                if "$ref" in ref
            ]
        else:
            ordered_group_keys = list(defs.keys())

        parts = []

        if section_title:
            parts.append(
                f'<h2 class="param-page-section">{html_module.escape(section_title)}</h2>\n'
            )

        for group_key in ordered_group_keys:
            group = defs.get(group_key)
            if group is None:
                continue

            properties = group.get("properties", {})
            if not properties:
                continue

            # Filter to visible or hidden parameters only
            visible_params = [
                (pname, pval)
                for pname, pval in properties.items()
                if bool(pval.get("hidden", False)) == show_hidden
            ]

            if not visible_params:
                continue

            required_params = set(group.get("required", []))
            title = group.get("title", group_key.replace("_", " ").title())
            description = group.get("description", "")
            fa_icon = group.get("fa_icon") or DEFAULT_ICON
            anchor = group_key.replace("_", "-")
            icon_html = f'<i class="{fa_icon}"></i>'

            parts.append(
                f'<div class="param-section" id="{anchor}">\n'
                f'<h3 class="param-section-header">{icon_html} {html_module.escape(title)}</h3>\n'
            )
            if description:
                parts.append(
                    f'<p class="param-section-desc">{html_module.escape(description)}</p>\n'
                )

            parts.append(
                '<table class="param-table">\n'
                "<thead><tr>"
                "<th>Parameter</th>"
                "<th>Default</th>"
                "<th>Description</th>"
                "</tr></thead>\n"
                "<tbody>\n"
            )

            for param_name, param in visible_params:
                is_required = param_name in required_params

                # --- Type badge ---
                raw_type = param.get("type", "")
                if isinstance(raw_type, list):
                    type_label = " | ".join(raw_type)
                else:
                    type_label = raw_type or "any"
                type_css = f"badge-type-{type_label.split()[0]}"
                type_badge = (
                    f'<span class="param-badge param-type {type_css}">'
                    f"{html_module.escape(type_label)}</span>"
                )

                # --- Badges under the parameter name ---
                badges = type_badge
                if is_required:
                    badges += ' <span class="param-badge param-required">required</span>'

                # --- Default value ---
                default = param.get("default")
                if default is None:
                    default_html = '<span class="param-no-default">-</span>'
                elif isinstance(default, bool):
                    val = "true" if default else "false"
                    default_html = f'<code class="param-default">{val}</code>'
                else:
                    default_html = (
                        f'<code class="param-default">'
                        f"{html_module.escape(str(default))}</code>"
                    )

                # --- Description ---
                desc_text = html_module.escape(param.get("description", ""))

                # --- Format badge (file-path, directory-path, etc.) ---
                fmt = param.get("format", "")
                fmt_html = ""
                if fmt:
                    fmt_html = (
                        f'<div class="param-format">'
                        f'<span class="param-badge param-format-badge">'
                        f"{html_module.escape(fmt)}</span></div>"
                    )

                # --- Enum chips ---
                enum_vals = param.get("enum", [])
                enum_html = ""
                if enum_vals:
                    chips = "".join(
                        f'<code class="param-enum-val">{html_module.escape(str(v))}</code> '
                        for v in enum_vals
                    )
                    enum_html = (
                        f'<div class="param-enum"><strong>Options:</strong> {chips}</div>'
                    )

                # --- Help text (collapsible <details>) ---
                help_text = param.get("help_text", "")
                help_html = ""
                if help_text:
                    help_html = (
                        '<details class="param-help">'
                        "<summary>More info</summary>"
                        f"<p>{html_module.escape(help_text)}</p>"
                        "</details>"
                    )

                parts.append(
                    f'<tr class="param-row">\n'
                    f'  <td class="param-name-cell">'
                    f'<code class="param-name">--{param_name}</code><br/>{badges}</td>\n'
                    f'  <td class="param-default-cell">{default_html}</td>\n'
                    f'  <td class="param-desc-cell">'
                    f"{desc_text}{fmt_html}{enum_html}{help_html}</td>\n"
                    f"</tr>\n"
                )

            parts.append("</tbody>\n</table>\n</div>\n\n")

        return "".join(parts)
