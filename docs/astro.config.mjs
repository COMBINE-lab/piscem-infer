import { defineConfig } from 'astro/config';
import starlight from '@astrojs/starlight';
import starlightThemeGalaxy from 'starlight-theme-galaxy';

export default defineConfig({
  site: 'https://combine-lab.github.io',
  base: '/piscem-infer',
  integrations: [
    starlight({
      plugins: [starlightThemeGalaxy()],
      title: 'piscem-infer',
      description:
        'Transcript-level abundance estimation from bulk-oriented RAD files, including multi-sample quantification with condition rescue.',
      social: [
        {
          icon: 'github',
          label: 'GitHub',
          href: 'https://github.com/COMBINE-lab/piscem-infer',
        },
      ],
      editLink: {
        baseUrl:
          'https://github.com/COMBINE-lab/piscem-infer/edit/multi-sample/docs/',
      },
      sidebar: [
        {
          label: 'Getting started',
          items: [
            { slug: 'getting-started/introduction' },
            { slug: 'getting-started/installation' },
            { slug: 'getting-started/quick-start' },
          ],
        },
        {
          label: 'Guides',
          items: [
            { slug: 'guides/building-inputs' },
            { slug: 'guides/single-sample-quantification' },
            { slug: 'guides/multi-sample-quantification' },
            { slug: 'guides/condition-rescue' },
          ],
        },
        {
          label: 'Reference',
          items: [
            { slug: 'reference/manifest' },
            { slug: 'reference/output-files' },
            { slug: 'reference/cli' },
          ],
        },
        {
          label: 'About',
          items: [{ slug: 'about/project-status' }],
        },
      ],
    }),
  ],
});
