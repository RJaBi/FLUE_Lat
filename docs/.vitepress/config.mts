// .vitepress/config.mts
import { withMermaid } from 'vitepress-plugin-mermaid'  // replace defineConfig import
import apiSidebar from '../api/_sidebar.json'


export default withMermaid({
  title: 'FLUE Documentation',
  base: '/FLUE_Lat/',
  head: [
  ['link', { rel: 'icon', href: '/images/lattice-063.png' }]
   ],
  markdown: {
    math: true,
    languages: ['fortran-free-form', 'fortran-fixed-form'],
    languageAlias: {
      'fortran': 'fortran-free-form',
      'f90': 'fortran-free-form',
      'f95': 'fortran-free-form',
      'f03': 'fortran-free-form',
      'f08': 'fortran-free-form',
      'f77': 'fortran-fixed-form',
      'F90': 'fortran-fixed-form',
    },
  },
  themeConfig: {
    nav: [
      { text: 'Home', link: '/' },
      { text: 'API', link: '/api/' },
    ],
    sidebar: {
      '/api/': [
        {
          text: 'API Reference',
          items: [
            { text: 'Overview', link: '/api/' },
          ],
        },
        ...apiSidebar,
      ],
    },
    search: {
      provider: 'local',
    },
  },
  mermaid: {},

})
